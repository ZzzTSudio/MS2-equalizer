#include "msequalizer.h"

#include "kissfft-master/kiss_fftr.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EQUALIZER_DEFAULT_RATE 8000
#define GAIN_ZERODB 1.0f

struct MSEqualizer {
    int rate;
    int nfft;
    float *frequency_response;
    float *fir;
    float *history;
    int fir_len;
    int history_len;
    int history_pos;
    int needs_update;
    int active;
};

static float clampf(float value, float min_value, float max_value) {
    if (value < min_value) return min_value;
    if (value > max_value) return max_value;
    return value;
}

static int equalizer_hz_to_index(const MSEqualizer *eq, int hz) {
    int ret;

    if (eq == NULL || hz < 0) return -1;
    if (hz > (eq->rate / 2)) hz = eq->rate / 2;

    ret = ((hz * eq->nfft) + (eq->rate / 2)) / eq->rate;
    if (ret == eq->nfft / 2) ret = (eq->nfft / 2) - 1;
    return ret;
}

static int equalizer_index_to_hz(const MSEqualizer *eq, int index) {
    return (index * eq->rate + eq->nfft / 2) / eq->nfft;
}

static void equalizer_flatten(MSEqualizer *eq) {
    size_t nfreqs;
    size_t i;

    nfreqs = (size_t)eq->nfft / 2;
    eq->frequency_response[0] = GAIN_ZERODB;
    for (i = 1; i < nfreqs; ++i) {
        eq->frequency_response[i] = GAIN_ZERODB;
    }
    eq->needs_update = 1;
}

static int equalizer_alloc_buffers(MSEqualizer *eq, int nfft) {
    float *frequency_response;
    float *fir;
    float *history;

    frequency_response = (float *)calloc((size_t)nfft / 2, sizeof(float));
    fir = (float *)calloc((size_t)nfft, sizeof(float));
    history = (float *)calloc((size_t)nfft > 0 ? (size_t)nfft - 1 : 0, sizeof(float));
    if (frequency_response == NULL || fir == NULL || history == NULL) {
        free(frequency_response);
        free(fir);
        free(history);
        return -1;
    }

    free(eq->frequency_response);
    free(eq->fir);
    free(eq->history);

    eq->frequency_response = frequency_response;
    eq->fir = fir;
    eq->history = history;
    eq->fir_len = nfft;
    eq->history_len = nfft - 1;
    eq->history_pos = 0;
    return 0;
}

static int equalizer_select_nfft(int rate) {
    if (rate < 16000) return 128;
    if (rate < 32000) return 256;
    return 512;
}

static int equalizer_rate_update(MSEqualizer *eq, int rate) {
    int nfft;

    if (eq == NULL || rate <= 0) return -1;

    nfft = equalizer_select_nfft(rate);
    if (equalizer_alloc_buffers(eq, nfft) != 0) return -1;

    eq->rate = rate;
    eq->nfft = nfft;
    equalizer_flatten(eq);
    return 0;
}

static float equalizer_compute_gainpoint(int f, int freq_0, float sqrt_gain, int freq_bw) {
    float k1;
    float k2;

    k1 = ((float)(f * f) - (float)(freq_0 * freq_0));
    k1 *= k1;
    k2 = (float)(f * freq_bw);
    k2 *= k2;
    return (k1 + k2 * sqrt_gain) / (k1 + k2 / sqrt_gain);
}

static void equalizer_point_set(MSEqualizer *eq, int index, float gain) {
    size_t nfreqs;

    nfreqs = (size_t)eq->nfft / 2;
    if (index >= 0 && (size_t)index < nfreqs) {
        eq->frequency_response[index] *= gain;
    }
}

static int equalizer_set_internal(MSEqualizer *eq, int freq_0, float gain, int freq_bw) {
    int i;
    int f;
    int delta_f;
    int mid;
    float current_gain;
    float sqrt_gain;

    if (eq == NULL || gain <= 0.0f) return -1;

    delta_f = equalizer_index_to_hz(eq, 1);
    sqrt_gain = sqrtf(gain);
    mid = equalizer_hz_to_index(eq, freq_0);
    if (mid < 0) return -1;

    freq_bw -= delta_f / 2;
    if (freq_bw < delta_f / 2) freq_bw = delta_f / 2;

    equalizer_point_set(eq, mid, gain);

    i = mid;
    current_gain = gain;
    do {
        ++i;
        f = equalizer_index_to_hz(eq, i);
        current_gain = equalizer_compute_gainpoint(f - delta_f, freq_0, sqrt_gain, freq_bw);
        equalizer_point_set(eq, i, current_gain);
    } while (i < eq->nfft / 2 && (current_gain > 1.1f || current_gain < 0.9f));

    i = mid;
    current_gain = gain;
    do {
        --i;
        f = equalizer_index_to_hz(eq, i);
        current_gain = equalizer_compute_gainpoint(f + delta_f, freq_0, sqrt_gain, freq_bw);
        equalizer_point_set(eq, i, current_gain);
    } while (i >= 0 && (current_gain > 1.1f || current_gain < 0.9f));

    eq->needs_update = 1;
    return 0;
}

static void time_shift(float *samples, int len) {
    int i;
    int half;

    half = len / 2;
    for (i = 0; i < half; ++i) {
        float tmp = samples[i];
        samples[i] = samples[i + half];
        samples[i + half] = tmp;
    }
}

static void norm_and_apodize(float *samples, int len) {
    int i;

    for (i = 0; i < len; ++i) {
        float x = (float)i * 2.0f * (float)M_PI / (float)len;
        float w = 0.54f - (0.46f * cosf(x));
        samples[i] *= w;
    }
}

static int equalizer_compute_impulse_response(MSEqualizer *eq) {
    kiss_fftr_cfg cfg;
    kiss_fft_cpx *freqdata;
    size_t bin_count;
    size_t i;

    cfg = kiss_fftr_alloc(eq->nfft, 1, NULL, NULL);
    if (cfg == NULL) return -1;

    bin_count = (size_t)eq->nfft / 2 + 1;
    freqdata = (kiss_fft_cpx *)calloc(bin_count, sizeof(kiss_fft_cpx));
    if (freqdata == NULL) {
        kiss_fftr_free(cfg);
        return -1;
    }

    freqdata[0].r = eq->frequency_response[0] / (float)eq->nfft;
    freqdata[0].i = 0.0f;
    for (i = 1; i < (size_t)eq->nfft / 2; ++i) {
        freqdata[i].r = eq->frequency_response[i] / (float)eq->nfft;
        freqdata[i].i = 0.0f;
    }
    freqdata[bin_count - 1].r = 1.0f / (float)eq->nfft;
    freqdata[bin_count - 1].i = 0.0f;

    kiss_fftri(cfg, freqdata, eq->fir);
    kiss_fftr_free(cfg);
    free(freqdata);

    time_shift(eq->fir, eq->fir_len);
    norm_and_apodize(eq->fir, eq->fir_len);
    memset(eq->history, 0, sizeof(float) * (size_t)eq->history_len);
    eq->history_pos = 0;
    eq->needs_update = 0;
    return 0;
}

static float equalizer_get_internal(const MSEqualizer *eq, int freqhz) {
    int idx;

    idx = equalizer_hz_to_index(eq, freqhz);
    if (idx < 0) return 0.0f;
    return eq->frequency_response[idx];
}

MSEqualizer *ms_equalizer_create(int sample_rate) {
    MSEqualizer *eq;

    if (sample_rate <= 0) sample_rate = EQUALIZER_DEFAULT_RATE;

    eq = (MSEqualizer *)calloc(1, sizeof(*eq));
    if (eq == NULL) return NULL;

    eq->active = 1;
    if (equalizer_rate_update(eq, sample_rate) != 0) {
        ms_equalizer_destroy(eq);
        return NULL;
    }

    return eq;
}

void ms_equalizer_destroy(MSEqualizer *eq) {
    if (eq == NULL) return;
    free(eq->frequency_response);
    free(eq->fir);
    free(eq->history);
    free(eq);
}

int ms_equalizer_set_sample_rate(MSEqualizer *eq, int sample_rate) {
    return equalizer_rate_update(eq, sample_rate);
}

int ms_equalizer_set_active(MSEqualizer *eq, int active) {
    if (eq == NULL) return -1;
    eq->active = active ? 1 : 0;
    return 0;
}

int ms_equalizer_is_active(const MSEqualizer *eq) {
    if (eq == NULL) return 0;
    return eq->active;
}

int ms_equalizer_set_gain(MSEqualizer *eq, const MSEqualizerGain *gain) {
    if (eq == NULL || gain == NULL) return -1;
    return equalizer_set_internal(eq, (int)gain->frequency, gain->gain, (int)gain->width);
}

int ms_equalizer_get_gain(const MSEqualizer *eq, MSEqualizerGain *gain) {
    if (eq == NULL || gain == NULL) return -1;
    gain->gain = equalizer_get_internal(eq, (int)gain->frequency);
    gain->width = 0.0f;
    return 0;
}

size_t ms_equalizer_get_num_frequencies(const MSEqualizer *eq) {
    if (eq == NULL) return 0;
    return (size_t)eq->nfft / 2;
}

int ms_equalizer_dump_state(const MSEqualizer *eq, float *table, size_t table_len) {
    size_t nfreqs;

    if (eq == NULL || table == NULL) return -1;

    nfreqs = ms_equalizer_get_num_frequencies(eq);
    if (table_len < nfreqs) return -1;

    memcpy(table, eq->frequency_response, sizeof(float) * nfreqs);
    return 0;
}

int ms_equalizer_process(MSEqualizer *eq, int16_t *samples, size_t nsamples) {
    size_t n;
    int taps;

    if (eq == NULL || samples == NULL) return -1;
    if (!eq->active || nsamples == 0) return 0;
    if (eq->needs_update && equalizer_compute_impulse_response(eq) != 0) return -1;

    taps = eq->fir_len;
    for (n = 0; n < nsamples; ++n) {
        float acc = eq->fir[0] * (float)samples[n];
        int hist_index = eq->history_pos;
        int tap;

        for (tap = 1; tap < taps; ++tap) {
            int idx = hist_index - (tap - 1);
            if (idx < 0) idx += eq->history_len;
            acc += eq->fir[tap] * eq->history[idx];
        }

        if (eq->history_len > 0) {
            eq->history_pos = (eq->history_pos + 1) % eq->history_len;
            eq->history[eq->history_pos] = (float)samples[n];
        }

        acc = clampf(acc, -32768.0f, 32767.0f);
        samples[n] = (int16_t)lrintf(acc);
    }

    return 0;
}

MSEqualizerGainArray ms_equalizer_parse_gains(const char *str) {
    MSEqualizerGainArray array = {0};

    if (str == NULL) return array;

    while (*str != '\0') {
        MSEqualizerGain gain;
        int bytes = 0;
        MSEqualizerGain *new_gains;

        if (sscanf(str, " %f:%f:%f %n", &gain.frequency, &gain.gain, &gain.width, &bytes) != 3) {
            break;
        }

        new_gains = (MSEqualizerGain *)realloc(array.gains, sizeof(*new_gains) * (array.count + 1));
        if (new_gains == NULL) {
            free(array.gains);
            array.gains = NULL;
            array.count = 0;
            return array;
        }

        array.gains = new_gains;
        array.gains[array.count++] = gain;
        str += bytes;
    }

    return array;
}

void ms_equalizer_free_gain_array(MSEqualizerGainArray *array) {
    if (array == NULL) return;
    free(array->gains);
    array->gains = NULL;
    array->count = 0;
}
