#include "msequalizer.h"

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void usage(const char *prog) {
    fprintf(stderr,
            "Usage: %s -r <sample_rate> [-b <samples_per_chunk>] [-g \"freq:gain:width ...\"] "
            "[-i <input.pcm>|-] [-o <output.pcm>|-]\n",
            prog);
}

static FILE *open_stream(const char *path, const char *mode, FILE *fallback) {
    if (path == NULL || strcmp(path, "-") == 0) return fallback;
    return fopen(path, mode);
}

int main(int argc, char **argv) {
    const char *input_path = "-";
    const char *output_path = "-";
    const char *gain_spec = NULL;
    int sample_rate = 16000;
    size_t chunk_samples = 320;
    int i;
    FILE *infile;
    FILE *outfile;
    int16_t *buffer;
    MSEqualizer *eq;
    MSEqualizerGainArray gains = {0};
    int exit_code = 1;

    for (i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            input_path = argv[++i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_path = argv[++i];
        } else if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
            sample_rate = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            chunk_samples = (size_t)strtoul(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            gain_spec = argv[++i];
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (sample_rate <= 0 || chunk_samples == 0) {
        usage(argv[0]);
        return 1;
    }

    infile = open_stream(input_path, "rb", stdin);
    if (infile == NULL) {
        fprintf(stderr, "failed to open input '%s': %s\n", input_path, strerror(errno));
        return 1;
    }

    outfile = open_stream(output_path, "wb", stdout);
    if (outfile == NULL) {
        fprintf(stderr, "failed to open output '%s': %s\n", output_path, strerror(errno));
        if (infile != stdin) fclose(infile);
        return 1;
    }

    eq = ms_equalizer_create(sample_rate);
    if (eq == NULL) {
        fprintf(stderr, "failed to create equalizer\n");
        goto cleanup_streams;
    }

    if (gain_spec != NULL) {
        size_t idx;

        gains = ms_equalizer_parse_gains(gain_spec);
        if (gains.count == 0) {
            fprintf(stderr, "failed to parse gain spec: %s\n", gain_spec);
            goto cleanup_equalizer;
        }

        for (idx = 0; idx < gains.count; ++idx) {
            if (ms_equalizer_set_gain(eq, &gains.gains[idx]) != 0) {
                fprintf(stderr, "failed to apply gain #%zu\n", idx);
                goto cleanup_equalizer;
            }
        }
    }

    buffer = (int16_t *)malloc(sizeof(*buffer) * chunk_samples);
    if (buffer == NULL) {
        fprintf(stderr, "failed to allocate audio buffer\n");
        goto cleanup_equalizer;
    }

    size_t read_count;

    while ((read_count = fread(buffer, sizeof(*buffer), chunk_samples, infile)) > 0) {
        if (ms_equalizer_process(eq, buffer, read_count) != 0) {
            fprintf(stderr, "equalizer processing failed\n");
            free(buffer);
            goto cleanup_equalizer;
        }
    
        if (fwrite(buffer, sizeof(*buffer), read_count, outfile) != read_count) {
            fprintf(stderr, "failed while writing output\n");
            free(buffer);
            goto cleanup_equalizer;
        }
    }
    
    if (ferror(infile)) {
        fprintf(stderr, "failed while reading input\n");
        free(buffer);
        goto cleanup_equalizer;
    }

    free(buffer);
    exit_code = 0;

cleanup_equalizer:
    ms_equalizer_free_gain_array(&gains);
    ms_equalizer_destroy(eq);
cleanup_streams:
    if (infile != stdin) fclose(infile);
    if (outfile != stdout) fclose(outfile);
    return exit_code;
}
