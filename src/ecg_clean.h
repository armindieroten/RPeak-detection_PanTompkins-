/*
 * ecg_clean.h
 *
 *  Created on: Jul 6, 2024
 *      Author: Hamidreza
 */

#ifndef ECG_CLEAN_H
#define ECG_CLEAN_H

#include <stdint.h>
#include <stddef.h>

void ecg_clean_nk(float* ecg_signal, size_t signal_length, uint32_t sampling_rate, float* clean_signal);
void ecg_clean_pt(float* ecg_signal, size_t signal_length, uint32_t sampling_rate, float* clean_signal);

#endif // ECG_CLEAN_H
