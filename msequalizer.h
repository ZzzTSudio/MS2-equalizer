#ifndef MSEQUALIZER_H
#define MSEQUALIZER_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MSEqualizer MSEqualizer;

typedef struct MSEqualizerGain {
    float frequency;
    float gain;
    float width;
} MSEqualizerGain;

typedef struct MSEqualizerGainArray {
    MSEqualizerGain *gains;
    size_t count;
} MSEqualizerGainArray;

MSEqualizer *ms_equalizer_create(int sample_rate);
void ms_equalizer_destroy(MSEqualizer *eq);

int ms_equalizer_set_sample_rate(MSEqualizer *eq, int sample_rate);
int ms_equalizer_set_active(MSEqualizer *eq, int active);
int ms_equalizer_is_active(const MSEqualizer *eq);

int ms_equalizer_set_gain(MSEqualizer *eq, const MSEqualizerGain *gain);
int ms_equalizer_get_gain(const MSEqualizer *eq, MSEqualizerGain *gain);

size_t ms_equalizer_get_num_frequencies(const MSEqualizer *eq);
int ms_equalizer_dump_state(const MSEqualizer *eq, float *table, size_t table_len);

int ms_equalizer_process(MSEqualizer *eq, int16_t *samples, size_t nsamples);

MSEqualizerGainArray ms_equalizer_parse_gains(const char *str);
void ms_equalizer_free_gain_array(MSEqualizerGainArray *array);

#ifdef __cplusplus
}
#endif

#endif
