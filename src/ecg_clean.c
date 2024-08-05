/*
 * ecg_clean.c
 *
 *  Created on: Jul 6, 2024
 *      Author: Hamidreza
 */


#include "ecg_clean.h"
#include <math.h>
#include <string.h>


// Define PI if not defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Function to compute Butterworth highpass filter coefficients
static void butterworth_highpass_coefficients(float* a, float* b, uint32_t sampling_rate, float cutoff, uint8_t order) {
  float nyquist = sampling_rate / 2.0f;
    float wc = tanf(M_PI * cutoff / sampling_rate);

    float a0 = 1.0f + sqrtf(2.0f) * wc + wc * wc;
    b[0] = (wc * wc) / a0;
    b[1] = -2.0f * (wc * wc) / a0;
    b[2] = (wc * wc) / a0;
    a[1] = 2.0f * (wc * wc - 1.0f) / a0;
    a[2] = (1.0f - sqrtf(2.0f) * wc + wc * wc) / a0;

    // Higher order coefficients (5th order)
    for (int i = 3; i <= 5; ++i) {
        b[i] = 0.0f;
        a[i] = 0.0f;
    }
}

// Function to apply Butterworth highpass filter
static void butterworth_highpass(const float* input_signal, float* output_signal, size_t signal_length, uint32_t sampling_rate, float cutoff, uint8_t order) {
     float a[6], b[6]; // Coefficient arrays for 5th order filter

    butterworth_highpass_coefficients(a, b, sampling_rate, cutoff,order);

    // Initialize history arrays for input and output
    float x[6] = {0.0f}; // Input history (x[n], x[n-1], ..., x[n-5])
    float y[6] = {0.0f}; // Output history (y[n], y[n-1], ..., y[n-5])

    // Apply filter to each sample
    for (size_t n = 0; n < signal_length; ++n) {
        // Update input history
        for (int i = 5; i > 0; --i) {
            x[i] = x[i - 1];
        }
        x[0] = input_signal[n];

        // Calculate output using difference equation
        float y0 = b[0] * x[0] + b[1] * x[1] + b[2] * x[2] - a[1] * y[1] - a[2] * y[2];

        // Update output history
        for (int i = 5; i > 0; --i) {
            y[i] = y[i - 1];
        }
        y[0] = y0;

        // Store filtered output
        output_signal[n] = y0;
    }
}



// Function to compute Butterworth bandpass filter coefficients
static void butterworth_bandpass_coefficients(float* a, float* b, uint32_t sampling_rate, float lowcut, float highcut, uint8_t order) {
    float nyquist = 0.5f * sampling_rate;
    float low = lowcut / nyquist;
    float high = highcut / nyquist;

    float omega_c = 2.0f * M_PI * ((low + high) / 2.0f);
    float bw = 2.0f * M_PI * (high - low);

    float alpha = sinf(bw / 2.0f) * sinh(log(2.0f) / 2.0f * bw / sinf(bw / 2.0f));
    float cos_omega_c = cosf(omega_c);

    // Calculate coefficients for a 5th order Butterworth bandpass filter
    float a0 = 1.0f + alpha;
    float a1 = -5.0f * cos_omega_c;
    float a2 = 10.0f;
    float a3 = -10.0f * cos_omega_c;
    float a4 = 5.0f;
    float a5 = -alpha;

    float b0 = alpha;
    float b1 = 0.0f;
    float b2 = -alpha;

    // Normalize coefficients
    b[0] = b0 / a0;
    b[1] = b1 / a0;
    b[2] = b2 / a0;
    b[3] = 0.0f; // 5th order filter has zeros for b3, b4, b5
    b[4] = 0.0f;
    b[5] = 0.0f;

    a[0] = 1.0f;
    a[1] = a1 / a0;
    a[2] = a2 / a0;
    a[3] = a3 / a0;
    a[4] = a4 / a0;
    a[5] = a5 / a0;
}

// Function to apply Butterworth bandpass filter
static void butterworth_bandpass(const float* input_signal, float* output_signal, size_t signal_length, uint32_t sampling_rate, float lowcut, float highcut, uint8_t order) {
    float a[6], b[6];
    butterworth_bandpass_coefficients(a, b, sampling_rate, lowcut, highcut, order);

    // Apply filter
    float x[6] = {0.0f};
    float y[6] = {0.0f};

    for (size_t i = 0; i < signal_length; ++i) {
        float x0 = input_signal[i];
        float y0 = b[0] * x0 + b[1] * x[1] + b[2] * x[2] + b[3] * x[3] + b[4] * x[4] + b[5] * x[5]
                   - a[1] * y[1] - a[2] * y[2] - a[3] * y[3] - a[4] * y[4] - a[5] * y[5];

        // Shift input and output arrays
        for (int j = order - 1; j > 0; --j) {
            x[j] = x[j - 1];
            y[j] = y[j - 1];
        }
        x[0] = x0;
        y[0] = y0;

        output_signal[i] = y0;
    }
}



// Function to apply powerline notch filter
static void powerline_filter(const float* input_signal, float* output_signal, size_t signal_length, uint32_t sampling_rate) {
    /** \brief
     *
     * \param
     * \param
     * \return
     *
     */
     /*
    float notch_freq = 50.0f; // Adjust this to 60.0f for regions with 60 Hz powerline frequency
    float w0 = 2.0f * M_PI * notch_freq / sampling_rate;
    float bw = w0 / 35.0f; // Bandwidth of the notch filter

    float a[3], b[3];
    b[0] = 1.0f;
    b[1] = -2.0f * cosf(w0);
    b[2] = 1.0f;
    a[0] = 1.0f + bw;
    a[1] = -2.0f * cosf(w0);
    a[2] = 1.0f - bw;

    // Normalize coefficients
    for (int i = 0; i < 3; ++i) {
        b[i] /= a[0];
        a[i] /= a[0];
    }

    // Apply filter
    float x1 = 0.0f, x2 = 0.0f, y1 = 0.0f, y2 = 0.0f;
    for (size_t i = 0; i < signal_length; ++i) {
        float x0 = input_signal[i];
        float y0 = b[0] * x0 + b[1] * x1 + b[2] * x2 - a[1] * y1 - a[2] * y2;

        output_signal[i] = y0;

        x2 = x1;
        x1 = x0;
        y2 = y1;
        y1 = y0;
    }
    */
    // FIR Notch filter coefficients (50 Hz notch)

}

// General signal filter function
void signal_filter(const float* input_signal, float* output_signal, size_t signal_length, uint32_t sampling_rate, float lowcut, float highcut,const char* method, uint8_t order) {
    if (strcmp(method, "butterworth") == 0) {
        if (lowcut > 0) {
            butterworth_highpass(input_signal, output_signal, signal_length, sampling_rate, lowcut,order);
        }
    }
    else if (strcmp(method, "butterworth-bp") == 0) {
            if (lowcut > 0 && highcut > lowcut) {
                butterworth_bandpass(input_signal, output_signal, signal_length, sampling_rate, lowcut,highcut, order);
            }
        }
    else if (strcmp(method, "powerline") == 0) {
        powerline_filter(input_signal, output_signal, signal_length, sampling_rate);
    }
}

// ECG signal cleaning function
void ecg_clean_nk(float* ecg_signal, size_t signal_length, uint32_t sampling_rate, float* clean_signal) {
    // Remove slow drift and dc offset with highpass Butterworth.
    signal_filter(ecg_signal, clean_signal, signal_length, sampling_rate, 0.5f, 0.0f,"butterworth", 5);

    // Remove powerline noise.
   signal_filter(clean_signal, clean_signal, signal_length, sampling_rate, 0.0f,0.0f, "powerline", 0);
}


// ECG signal cleaning function
void ecg_clean_pt(float* ecg_signal, size_t signal_length, uint32_t sampling_rate, float* clean_signal) {
    // Remove slow drift and dc offset with highpass Butterworth.
    signal_filter(ecg_signal, clean_signal, signal_length, sampling_rate, 5.0f, 15.0f,"butterworth-bp", 5);

    // Remove powerline noise.
   // signal_filter(clean_signal, clean_signal, signal_length, sampling_rate, 0.0f,0.0f, "powerline", 0);
}
