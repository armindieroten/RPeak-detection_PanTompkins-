
#include <stdio.h>
#include <stdlib.h>


#include "ecg_clean.h"
#include "filter.h"
#include "PanTompkins.h"

//static uint8_t data_array[4000] = {184, 30, 133, 190, 184, 30, 133, 190, 184, 30, 133, 190, 184, 30, 133, 190, 184, 30, 133, 190, 184, 30, 133, 190, 184, 30, 133, 190, 184, 30, 133, 190, 184, 30, 133, 190, 0, 0, 128, 190, 143, 194, 117, 190, 102, 102, 102, 190, 61, 10, 87, 190, 205, 204, 76, 190, 205, 204, 76, 190, 205, 204, 76, 190, 236, 81, 56, 190, 123, 20, 46, 190, 82, 184, 30, 190, 195, 245, 40, 190, 123, 20, 46, 190, 51, 51, 51, 190, 123, 20, 46, 190, 195, 245, 40, 190, 10, 215, 35, 190, 123, 20, 46, 190, 51, 51, 51, 190, 164, 112, 61, 190, 51, 51, 51, 190, 123, 20, 46, 190, 195, 245, 40, 190, 51, 51, 51, 190, 205, 204, 76, 190, 133, 235, 81, 190, 61, 10, 87, 190, 61, 10, 87, 190, 61, 10, 87, 190, 102, 102, 102, 190, 143, 194, 117, 190, 0, 0, 128, 190, 72, 225, 122, 190, 72, 225, 122, 190, 0, 0, 128, 190, 92, 143, 130, 190, 113, 61, 138, 190, 225, 122, 148, 190, 246, 40, 156, 190, 82, 184, 158, 190, 246, 40, 156, 190, 246, 40, 156, 190, 82, 184, 158, 190, 82, 184, 158, 190, 82, 184, 158, 190, 61, 10, 151, 190, 154, 153, 153, 190, 10, 215, 163, 190, 102, 102, 166, 190, 102, 102, 166, 190, 10, 215, 163, 190, 246, 40, 156, 190, 225, 122, 148, 190, 61, 10, 151, 190, 225, 122, 148, 190, 246, 40, 156, 190, 61, 10, 151, 190, 133, 235, 145, 190, 61, 10, 151, 190, 61, 10, 151, 190, 154, 153, 153, 190, 246, 40, 156, 190, 133, 235, 145, 190, 133, 235, 145, 190, 41, 92, 143, 190, 133, 235, 145, 190, 225, 122, 148, 190, 61, 10, 151, 190, 61, 10, 151, 190, 41, 92, 143, 190, 205, 204, 140, 190, 20, 174, 135, 190, 92, 143, 130, 190, 184, 30, 133, 190, 72, 225, 122, 190, 215, 163, 112, 190, 215, 163, 112, 190, 143, 194, 117, 190, 143, 194, 117, 190, 92, 143, 130, 190, 0, 0, 128, 190, 72, 225, 122, 190, 72, 225, 122, 190, 92, 143, 130, 190, 184, 30, 133, 190, 20, 174, 135, 190, 184, 30, 133, 190, 113, 61, 138, 190, 205, 204, 140, 190, 225, 122, 148, 190, 154, 153, 153, 190, 246, 40, 156, 190, 82, 184, 158, 190, 246, 40, 156, 190, 225, 122, 148, 190, 133, 235, 145, 190, 154, 153, 153, 190, 154, 153, 153, 190, 61, 10, 151, 190, 133, 235, 145, 190, 133, 235, 145, 190, 41, 92, 143, 190, 225, 122, 148, 190, 246, 40, 156, 190, 154, 153, 153, 190, 154, 153, 153, 190, 246, 40, 156, 190, 174, 71, 161, 190, 31, 133, 171, 190, 215, 163, 176, 190, 123, 20, 174, 190, 195, 245, 168, 190, 102, 102, 166, 190, 102, 102, 166, 190, 195, 245, 168, 190, 102, 102, 166, 190, 174, 71, 161, 190, 154, 153, 153, 190, 246, 40, 156, 190, 174, 71, 161, 190, 31, 133, 171, 190, 31, 133, 171, 190, 102, 102, 166, 190, 82, 184, 158, 190, 246, 40, 156, 190, 10, 215, 163, 190, 195, 245, 168, 190, 195, 245, 168, 190, 82, 184, 158, 190, 246, 40, 156, 190, 82, 184, 158, 190, 10, 215, 163, 190, 31, 133, 171, 190, 195, 245, 168, 190, 174, 71, 161, 190, 246, 40, 156, 190, 20, 174, 135, 190, 143, 194, 117, 190, 61, 10, 87, 190, 92, 143, 66, 190, 123, 20, 46, 190, 184, 30, 5, 190, 10, 215, 163, 189, 10, 215, 163, 59, 205, 204, 204, 61, 164, 112, 61, 62, 246, 40, 156, 62, 246, 40, 220, 62, 225, 122, 20, 63, 225, 122, 52, 63, 10, 215, 67, 63, 143, 194, 85, 63, 10, 215, 99, 63, 82, 184, 94, 63, 0, 0, 64, 63, 51, 51, 19, 63, 0, 0, 192, 62, 154, 153, 25, 62, 143, 194, 245, 188, 184, 30, 5, 190, 236, 81, 56, 190, 174, 71, 97, 190, 31, 133, 107, 190, 215, 163, 112, 190, 246, 40, 92, 190, 246, 40, 92, 190, 174, 71, 97, 190, 215, 163, 112, 190, 92, 143, 130, 190, 20, 174, 135, 190, 205, 204, 140, 190, 205, 204, 140, 190, 61, 10, 151, 190, 154, 153, 153, 190, 246, 40, 156, 190, 82, 184, 158, 190, 154, 153, 153, 190, 174, 71, 161, 190, 102, 102, 166, 190, 215, 163, 176, 190, 31, 133, 171, 190, 123, 20, 174, 190, 10, 215, 163, 190, 61, 10, 151, 190, 133, 235, 145, 190, 61, 10, 151, 190, 174, 71, 161, 190, 10, 215, 163, 190, 102, 102, 166, 190, 102, 102, 166, 190, 102, 102, 166, 190, 31, 133, 171, 190, 102, 102, 166, 190, 31, 133, 171, 190, 10, 215, 163, 190, 10, 215, 163, 190, 82, 184, 158, 190, 82, 184, 158, 190, 154, 153, 153, 190, 41, 92, 143, 190, 92, 143, 130, 190, 215, 163, 112, 190, 61, 10, 87, 190, 133, 235, 81, 190, 205, 204, 76, 190, 133, 235, 81, 190, 246, 40, 92, 190, 20, 174, 71, 190, 20, 174, 71, 190, 133, 235, 81, 190, 61, 10, 87, 190, 246, 40, 92, 190, 61, 10, 87, 190, 20, 174, 71, 190, 51, 51, 51, 190, 195, 245, 40, 190, 51, 51, 51, 190, 195, 245, 40, 190, 10, 215, 35, 190, 10, 215, 35, 190, 10, 215, 35, 190, 195, 245, 40, 190, 123, 20, 46, 190, 164, 112, 61, 190, 205, 204, 76, 190, 92, 143, 66, 190, 236, 81, 56, 190, 51, 51, 51, 190, 195, 245, 40, 190, 51, 51, 51, 190, 123, 20, 46, 190, 195, 245, 40, 190, 123, 20, 46, 190, 51, 51, 51, 190, 92, 143, 66, 190, 20, 174, 71, 190, 133, 235, 81, 190, 20, 174, 71, 190, 61, 10, 87, 190, 133, 235, 81, 190, 246, 40, 92, 190, 102, 102, 102, 190, 102, 102, 102, 190, 174, 71, 97, 190, 31, 133, 107, 190, 0, 0, 128, 190, 20, 174, 135, 190, 225, 122, 148, 190, 133, 235, 145, 190, 205, 204, 140, 190, 205, 204, 140, 190, 154, 153, 153, 190, 195, 245, 168, 190, 31, 133, 171, 190, 82, 184, 158, 190, 41, 92, 143, 190, 20, 174, 135, 190, 92, 143, 130, 190, 0, 0, 128, 190, 92, 143, 130, 190, 113, 61, 138, 190, 184, 30, 133, 190, 0, 0, 128, 190, 143, 194, 117, 190, 143, 194, 117, 190, 143, 194, 117, 190, 246, 40, 92, 190, 246, 40, 92, 190, 174, 71, 97, 190, 215, 163, 112, 190, 143, 194, 117, 190, 184, 30, 133, 190, 92, 143, 130, 190, 92, 143, 130, 190, 184, 30, 133, 190, 92, 143, 130, 190, 184, 30, 133, 190, 20, 174, 135, 190, 92, 143, 130, 190, 92, 143, 130, 190, 20, 174, 135, 190, 113, 61, 138, 190, 205, 204, 140, 190, 41, 92, 143, 190, 20, 174, 135, 190, 113, 61, 138, 190, 20, 174, 135, 190, 205, 204, 140, 190, 205, 204, 140, 190, 225, 122, 148, 190, 61, 10, 151, 190, 154, 153, 153, 190, 225, 122, 148, 190, 61, 10, 151, 190, 246, 40, 156, 190, 174, 71, 161, 190, 174, 71, 161, 190, 82, 184, 158, 190, 246, 40, 156, 190, 246, 40, 156, 190, 61, 10, 151, 190, 10, 215, 163, 190, 195, 245, 168, 190, 195, 245, 168, 190, 195, 245, 168, 190, 31, 133, 171, 190, 51, 51, 179, 190, 123, 20, 174, 190, 102, 102, 166, 190, 246, 40, 156, 190, 61, 10, 151, 190, 154, 153, 153, 190, 10, 215, 163, 190, 51, 51, 179, 190, 236, 81, 184, 190, 215, 163, 176, 190, 174, 71, 161, 190, 225, 122, 148, 190, 225, 122, 148, 190, 61, 10, 151, 190, 31, 133, 171, 190, 31, 133, 171, 190, 10, 215, 163, 190, 61, 10, 151, 190, 225, 122, 148, 190, 133, 235, 145, 190, 205, 204, 140, 190, 41, 92, 143, 190, 113, 61, 138, 190, 20, 174, 135, 190, 133, 235, 145, 190, 133, 235, 145, 190, 113, 61, 138, 190, 205, 204, 140, 190, 225, 122, 148, 190, 154, 153, 153, 190, 82, 184, 158, 190, 174, 71, 161, 190, 10, 215, 163, 190, 246, 40, 156, 190, 225, 122, 148, 190, 41, 92, 143, 190, 133, 235, 145, 190, 225, 122, 148, 190, 225, 122, 148, 190, 225, 122, 148, 190, 61, 10, 151, 190, 246, 40, 156, 190, 82, 184, 158, 190, 225, 122, 148, 190, 113, 61, 138, 190, 92, 143, 130, 190, 0, 0, 128, 190, 41, 92, 143, 190, 174, 71, 161, 190, 51, 51, 179, 190, 72, 225, 186, 190, 236, 81, 184, 190, 215, 163, 176, 190, 31, 133, 171, 190, 215, 163, 176, 190, 102, 102, 166, 190, 82, 184, 158, 190, 61, 10, 151, 190, 154, 153, 153, 190, 133, 235, 145, 190, 113, 61, 138, 190, 20, 174, 135, 190, 184, 30, 133, 190, 225, 122, 148, 190, 225, 122, 148, 190, 133, 235, 145, 190, 20, 174, 135, 190, 20, 174, 135, 190, 143, 194, 117, 190, 31, 133, 107, 190, 143, 194, 117, 190, 113, 61, 138, 190, 41, 92, 143, 190, 133, 235, 145, 190, 113, 61, 138, 190, 72, 225, 122, 190, 102, 102, 102, 190, 246, 40, 92, 190, 246, 40, 92, 190, 215, 163, 112, 190, 72, 225, 122, 190, 92, 143, 130, 190, 184, 30, 133, 190, 205, 204, 140, 190, 61, 10, 151, 190, 154, 153, 153, 190, 41, 92, 143, 190, 20, 174, 135, 190, 215, 163, 112, 190, 61, 10, 87, 190, 20, 174, 71, 190, 236, 81, 56, 190, 143, 194, 245, 189, 10, 215, 35, 189, 205, 204, 204, 61, 215, 163, 112, 62, 51, 51, 179, 62, 51, 51, 243, 62, 82, 184, 30, 63, 246, 40, 60, 63, 154, 153, 89, 63, 133, 235, 113, 63, 72, 225, 122, 63, 10, 215, 99, 63, 236, 81, 56, 63, 0, 0, 0, 63, 113, 61, 138, 62, 236, 81, 56, 61, 143, 194, 245, 189, 72, 225, 122, 190, 195, 245, 168, 190, 51, 51, 179, 190, 154, 153, 153, 190, 215, 163, 112, 190, 133, 235, 81, 190, 20, 174, 71, 190, 205, 204, 76, 190, 20, 174, 71, 190, 205, 204, 76, 190, 61, 10, 87, 190, 246, 40, 92, 190, 31, 133, 107, 190, 215, 163, 112, 190, 31, 133, 107, 190, 174, 71, 97, 190, 102, 102, 102, 190, 143, 194, 117, 190, 0, 0, 128, 190, 215, 163, 112, 190, 143, 194, 117, 190, 72, 225, 122, 190, 72, 225, 122, 190, 92, 143, 130, 190, 92, 143, 130, 190, 92, 143, 130, 190, 0, 0, 128, 190, 0, 0, 128, 190, 0, 0, 128, 190, 143, 194, 117, 190, 143, 194, 117, 190, 215, 163, 112, 190, 174, 71, 97, 190, 102, 102, 102, 190, 31, 133, 107, 190, 174, 71, 97, 190, 133, 235, 81, 190, 92, 143, 66, 190, 236, 81, 56, 190, 236, 81, 56, 190, 164, 112, 61, 190, 246, 40, 92, 190, 61, 10, 87, 190, 236, 81, 56, 190, 195, 245, 40, 190, 41, 92, 15, 190, 0, 0, 0, 190, 184, 30, 5, 190, 143, 194, 245, 189, 143, 194, 245, 189, 174, 71, 225, 189, 174, 71, 225, 189, 61, 10, 215, 189, 205, 204, 204, 189, 92, 143, 194, 189, 236, 81, 184, 189, 123, 20, 174, 189, 154, 153, 153, 189, 123, 20, 174, 189, 123, 20, 174, 189, 236, 81, 184, 189, 184, 30, 133, 189, 205, 204, 76, 189, 236, 81, 56, 189, 174, 71, 97, 189, 143, 194, 117, 189, 143, 194, 117, 189, 236, 81, 56, 189, 236, 81, 56, 189, 205, 204, 204, 188, 205, 204, 204, 188, 143, 194, 245, 188, 205, 204, 204, 188, 205, 204, 204, 188, 10, 215, 35, 188, 10, 215, 163, 187, 0, 0, 0, 0, 10, 215, 35, 188, 10, 215, 163, 187, 10, 215, 163, 187, 10, 215, 35, 188, 143, 194, 117, 188, 143, 194, 245, 188, 205, 204, 76, 189, 154, 153, 153, 189, 143, 194, 117, 189, 41, 92, 143, 189, 41, 92, 143, 189, 10, 215, 163, 189, 236, 81, 184, 189, 174, 71, 225, 189, 236, 81, 184, 189, 92, 143, 194, 189, 92, 143, 194, 189, 174, 71, 225, 189, 41, 92, 15, 190, 225, 122, 20, 190, 10, 215, 35, 190, 225, 122, 20, 190, 82, 184, 30, 190, 154, 153, 25, 190, 10, 215, 35, 190, 10, 215, 35, 190, 10, 215, 35, 190, 154, 153, 25, 190, 41, 92, 15, 190, 154, 153, 25, 190, 82, 184, 30, 190, 123, 20, 46, 190, 154, 153, 25, 190, 41, 92, 15, 190, 0, 0, 0, 190, 184, 30, 5, 190, 41, 92, 15, 190, 82, 184, 30, 190, 41, 92, 15, 190, 0, 0, 0, 190, 184, 30, 5, 190, 113, 61, 10, 190, 225, 122, 20, 190, 82, 184, 30, 190, 225, 122, 20, 190, 184, 30, 5, 190, 0, 0, 0, 190, 113, 61, 10, 190, 113, 61, 10, 190, 113, 61, 10, 190, 0, 0, 0, 190, 174, 71, 225, 189, 174, 71, 225, 189, 92, 143, 194, 189, 31, 133, 235, 189, 31, 133, 235, 189, 174, 71, 225, 189, 92, 143, 194, 189, 61, 10, 215, 189, 31, 133, 235, 189, 184, 30, 5, 190, 113, 61, 10, 190, 0, 0, 0, 190, 184, 30, 5, 190, 0, 0, 0, 190, 41, 92, 15, 190, 154, 153, 25, 190, 154, 153, 25, 190, 225, 122, 20, 190, 113, 61, 10, 190, 113, 61, 10, 190, 82, 184, 30, 190, 10, 215, 35, 190, 51, 51, 51, 190, 82, 184, 30, 190, 154, 153, 25, 190, 82, 184, 30, 190, 10, 215, 35, 190, 51, 51, 51, 190, 51, 51, 51, 190, 10, 215, 35, 190, 10, 215, 35, 190, 10, 215, 35, 190, 10, 215, 35, 190, 10, 215, 35, 190, 123, 20, 46, 190, 82, 184, 30, 190, 225, 122, 20, 190, 82, 184, 30, 190, 51, 51, 51, 190, 123, 20, 46, 190, 123, 20, 46, 190, 82, 184, 30, 190, 41, 92, 15, 190, 225, 122, 20, 190, 10, 215, 35, 190, 236, 81, 56, 190, 236, 81, 56, 190, 51, 51, 51, 190, 82, 184, 30, 190, 82, 184, 30, 190, 82, 184, 30, 190, 82, 184, 30, 190, 195, 245, 40, 190, 10, 215, 35, 190, 225, 122, 20, 190, 225, 122, 20, 190, 82, 184, 30, 190, 154, 153, 25, 190, 82, 184, 30, 190, 41, 92, 15, 190, 113, 61, 10, 190, 113, 61, 10, 190, 154, 153, 25, 190, 10, 215, 35, 190, 51, 51, 51, 190, 195, 245, 40, 190, 225, 122, 20, 190, 41, 92, 15, 190, 82, 184, 30, 190, 51, 51, 51, 190, 92, 143, 66, 190, 51, 51, 51, 190, 195, 245, 40, 190, 10, 215, 35, 190, 236, 81, 56, 190, 20, 174, 71, 190, 205, 204, 76, 190, 164, 112, 61, 190, 51, 51, 51, 190, 236, 81, 56, 190, 164, 112, 61, 190, 92, 143, 66, 190, 92, 143, 66, 190, 236, 81, 56, 190, 236, 81, 56, 190, 123, 20, 46, 190, 195, 245, 40, 190, 236, 81, 56, 190, 92, 143, 66, 190, 123, 20, 46, 190, 10, 215, 35, 190, 195, 245, 40, 190, 51, 51, 51, 190, 51, 51, 51, 190, 51, 51, 51, 190, 123, 20, 46, 190, 123, 20, 46, 190, 123, 20, 46, 190, 51, 51, 51, 190, 236, 81, 56, 190, 92, 143, 66, 190, 123, 20, 46, 190, 10, 215, 35, 190, 82, 184, 30, 190, 10, 215, 35, 190, 123, 20, 46, 190, 195, 245, 40, 190, 154, 153, 25, 190, 113, 61, 10, 190, 225, 122, 20, 190, 123, 20, 46, 190, 164, 112, 61, 190, 164, 112, 61, 190, 195, 245, 40, 190, 123, 20, 46, 190, 10, 215, 35, 190, 123, 20, 46, 190, 236, 81, 56, 190, 20, 174, 71, 190, 236, 81, 56, 190, 41, 92, 15, 190, 0, 0, 0, 190, 61, 10, 215, 189, 92, 143, 194, 189, 154, 153, 153, 189, 41, 92, 15, 189, 205, 204, 204, 60, 31, 133, 235, 61, 174, 71, 97, 62, 102, 102, 166, 62, 174, 71, 225, 62, 143, 194, 21, 63, 143, 194, 53, 63, 133, 235, 81, 63, 113, 61, 106, 63, 215, 163, 128, 63, 61, 10, 135, 63, 184, 30, 133, 63, 41, 92, 111, 63, 92, 143, 66, 63, 123, 20, 14, 63, 195, 245, 168, 62, 113, 61, 10, 62, 0, 0, 0, 0, 92, 143, 194, 189, 82, 184, 30, 190, 20, 174, 71, 190, 92, 143, 66, 190, 82, 184, 30, 190, 31, 133, 235, 189, 61, 10, 215, 189, 61, 10, 215, 189, 184, 30, 5, 190, 225, 122, 20, 190, 195, 245, 40, 190, 195, 245, 40, 190, 10, 215, 35, 190, 10, 215, 35, 190, 123, 20, 46, 190, 92, 143, 66, 190, 92, 143, 66, 190, 236, 81, 56, 190, 123, 20, 46, 190, 51, 51, 51, 190, 164, 112, 61, 190, 236, 81, 56, 190, 92, 143, 66, 190, 164, 112, 61, 190, 51, 51, 51, 190, 51, 51, 51, 190, 123, 20, 46, 190, 51, 51, 51, 190, 164, 112, 61, 190, 236, 81, 56, 190, 82, 184, 30, 190, 154, 153, 25, 190, 10, 215, 35, 190, 195, 245, 40, 190, 123, 20, 46, 190, 10, 215, 35, 190, 225, 122, 20, 190, 82, 184, 30, 190, 10, 215, 35, 190, 123, 20, 46, 190, 195, 245, 40, 190, 10, 215, 35, 190, 154, 153, 25, 190, 225, 122, 20, 190, 225, 122, 20, 190, 154, 153, 25, 190, 82, 184, 30, 190, 225, 122, 20, 190, 41, 92, 15, 190, 113, 61, 10, 190, 225, 122, 20, 190, 225, 122, 20, 190, 154, 153, 25, 190, 113, 61, 10, 190, 143, 194, 245, 189, 0, 0, 0, 190, 0, 0, 0, 190, 113, 61, 10, 190, 41, 92, 15, 190, 0, 0, 0, 190, 31, 133, 235, 189, 184, 30, 5, 190, 0, 0, 0, 190, 113, 61, 10, 190, 184, 30, 5, 190, 0, 0, 0, 190, 61, 10, 215, 189, 205, 204, 204, 189, 205, 204, 204, 189, 236, 81, 184, 189, 10, 215, 163, 189, 174, 71, 97, 189, 10, 215, 35, 189, 143, 194, 245, 188, 10, 215, 35, 189, 236, 81, 56, 189, 10, 215, 35, 189, 205, 204, 204, 188, 10, 215, 163, 188, 10, 215, 35, 188, 143, 194, 117, 188, 205, 204, 204, 188, 41, 92, 15, 189, 205, 204, 204, 188, 143, 194, 245, 188, 10, 215, 35, 189, 174, 71, 97, 189, 154, 153, 153, 189, 123, 20, 174, 189, 123, 20, 174, 189, 154, 153, 153, 189, 236, 81, 184, 189, 205, 204, 204, 189, 184, 30, 5, 190, 113, 61, 10, 190, 41, 92, 15, 190, 41, 92, 15, 190, 225, 122, 20, 190, 154, 153, 25, 190, 195, 245, 40, 190, 123, 20, 46, 190, 195, 245, 40, 190, 82, 184, 30, 190, 123, 20, 46, 190, 123, 20, 46, 190, 51, 51, 51, 190, 51, 51, 51, 190, 123, 20, 46, 190, 82, 184, 30, 190, 10, 215, 35, 190, 195, 245, 40, 190, 51, 51, 51, 190, 123, 20, 46, 190, 10, 215, 35, 190, 10, 215, 35, 190, 123, 20, 46, 190, 123, 20, 46, 190, 92, 143, 66, 190, 164, 112, 61, 190, 92, 143, 66, 190, 164, 112, 61, 190, 92, 143, 66, 190, 164, 112, 61, 190, 20, 174, 71, 190, 205, 204, 76, 190, 20, 174, 71, 190, 20, 174, 71, 190, 164, 112, 61, 190, 236, 81, 56, 190, 236, 81, 56, 190, 195, 245, 40, 190, 123, 20, 46, 190, 10, 215, 35, 190, 195, 245, 40, 190, 10, 215, 35, 190, 10, 215, 35, 190, 10, 215, 35, 190, 154, 153, 25, 190, 41, 92, 15, 190, 225, 122, 20, 190, 225, 122, 20, 190, 10, 215, 35, 190, 51, 51, 51, 190, 195, 245, 40, 190, 10, 215, 35, 190, 195, 245, 40, 190, 51, 51, 51, 190, 123, 20, 46, 190, 164, 112, 61, 190, 51, 51, 51, 190, 123, 20, 46, 190, 51, 51, 51, 190, 236, 81, 56, 190, 20, 174, 71, 190, 20, 174, 71, 190, 164, 112, 61, 190, 236, 81, 56, 190, 195, 245, 40, 190, 236, 81, 56, 190, 92, 143, 66, 190, 20, 174, 71, 190, 20, 174, 71, 190, 164, 112, 61, 190, 236, 81, 56, 190, 92, 143, 66, 190, 20, 174, 71, 190, 61, 10, 87, 190, 20, 174, 71, 190, 92, 143, 66, 190, 92, 143, 66, 190, 133, 235, 81, 190, 92, 143, 66, 190, 61, 10, 87, 190, 92, 143, 66, 190, 92, 143, 66, 190, 20, 174, 71, 190, 133, 235, 81, 190, 133, 235, 81, 190, 246, 40, 92, 190, 61, 10, 87, 190, 20, 174, 71, 190, 20, 174, 71, 190, 133, 235, 81, 190, 133, 235, 81, 190, 246, 40, 92, 190, 174, 71, 97, 190, 61, 10, 87, 190, 174, 71, 97, 190, 174, 71, 97, 190, 31, 133, 107, 190, 174, 71, 97, 190, 205, 204, 76, 190, 164, 112, 61, 190, 195, 245, 40, 190, 154, 153, 25, 190, 41, 92, 15, 190, 143, 194, 245, 189, 205, 204, 204, 189, 123, 20, 174, 189, 143, 194, 117, 189, 0, 0, 0, 0, 154, 153, 153, 61, 195, 245, 40, 62, 113, 61, 138, 62, 72, 225, 186, 62, 51, 51, 243, 62, 154, 153, 25, 63, 215, 163, 48, 63, 184, 30, 69, 63, 154, 153, 89, 63, 195, 245, 104, 63, 195, 245, 104, 63, 143, 194, 85, 63, 154, 153, 57, 63, 61, 10, 23, 63, 246, 40, 220, 62, 0, 0, 128, 62, 205, 204, 204, 61, 0, 0, 0, 0, 184, 30, 133, 189, 61, 10, 215, 189, 143, 194, 245, 189, 31, 133, 235, 189, 61, 10, 215, 189, 174, 71, 225, 189, 113, 61, 10, 190, 82, 184, 30, 190, 195, 245, 40, 190, 123, 20, 46, 190, 195, 245, 40, 190, 236, 81, 56, 190, 20, 174, 71, 190, 205, 204, 76, 190, 205, 204, 76, 190, 205, 204, 76, 190, 61, 10, 87, 190, 246, 40, 92, 190, 61, 10, 87, 190, 246, 40, 92, 190, 246, 40, 92, 190, 61, 10, 87, 190, 102, 102, 102, 190, 102, 102, 102, 190, 143, 194, 117, 190, 143, 194, 117, 190, 72, 225, 122, 190, 143, 194, 117, 190, 143, 194, 117, 190, 0, 0, 128, 190, 184, 30, 133, 190, 184, 30, 133, 190, 0, 0, 128, 190, 92, 143, 130, 190, 92, 143, 130, 190, 20, 174, 135, 190, 205, 204, 140, 190, 225, 122, 148, 190, 133, 235, 145, 190, 133, 235, 145, 190, 61, 10, 151, 190, 82, 184, 158, 190, 174, 71, 161, 190, 154, 153, 153, 190, 20, 174, 135, 190, 0, 0, 128, 190, 215, 163, 112, 190, 215, 163, 112, 190, 215, 163, 112, 190, 102, 102, 102, 190, 133, 235, 81, 190, 236, 81, 56, 190, 123, 20, 46, 190, 123, 20, 46, 190, 51, 51, 51, 190, 236, 81, 56, 190, 195, 245, 40, 190, 225, 122, 20, 190, 154, 153, 25, 190, 225, 122, 20, 190, 113, 61, 10, 190, 41, 92, 15, 190, 41, 92, 15, 190, 113, 61, 10, 190, 0, 0, 0, 190, 0, 0, 0, 190, 143, 194, 245, 189, 31, 133, 235, 189, 31, 133, 235, 189, 61, 10, 215, 189, 205, 204, 204, 189, 61, 10, 215, 189, 61, 10, 215, 189, 174, 71, 225, 189, 205, 204, 204, 189, 92, 143, 194, 189, 123, 20, 174, 189, 205, 204, 204, 189, 205, 204, 204, 189, 61, 10, 215, 189, 92, 143, 194, 189, 236, 81, 184, 189, 92, 143, 194, 189, 61, 10, 215, 189, 31, 133, 235, 189, 0, 0, 0, 190};


#define SIGNAL_SIZE 1000
#define SAMPLING_RATE 360
//#define PEAKS_SIZE 100
//#define BUFFER_SIZE (SIGNAL_SIZE + PEAKS_SIZE)

#define BUFFER_SIZE (SIGNAL_SIZE)

uint8_t rx_buffer[BUFFER_SIZE];
float *ecg_signal;
//float r_peaks[PEAKS_SIZE];

//float ecg_signal[10000];
float clean_signal[SIGNAL_SIZE];
int ecg_peaks[SIGNAL_SIZE];


void write_array_to_file(const char *filename, float *array, size_t length) {
    FILE *file = fopen(filename, "wb");
    if (file == NULL) {
        perror("Unable to open file");
        return;
    }

    fwrite(array, sizeof(float), length, file);
    fclose(file);
}


int append_to_file(const char *filename, int value) {
    FILE *file;

    // Open the file in append mode
    file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }

    // Write the integer to the file
    fprintf(file, "%d\n", value);

    // Close the file
    fclose(file);

    return 0;
}

int read_floats_from_file(const char *filename, float **array, int read_count);

int main()
{
    uint32_t read_count = 500000;
    const char *filename = "ecg.txt";

    if (read_floats_from_file(filename, &ecg_signal, read_count) != 0) {
        fprintf(stderr, "Error reading floats from file\n");
        return EXIT_FAILURE;
    }
    // Free the allocated memory

    //float *ecg_signal = (float *)data_array;

    uint32_t beat_mask[500000]={0};

    int16_t delay, Rcount,s1, s2, s3, s4, s5, ThF1;
	uint16_t ThI1, SPKI, NPKI;
	int32_t RLoc, c, SampleCount;
	SampleCount = 0;

	Rcount = 0;
	errno_t err, err1;

    PT_init();
    // ------ Pass the signal sample by sample mimicing a real-time scenario ----------- //
	for (uint32_t i = 0; i < read_count; i++) {
		++SampleCount;
		c = ecg_signal[i]*100;

		delay = PT_StateMachine((int16_t) c);							// This is the main function of the algorithm

		// ------- A positive delay to current sample is returned in case of beat detection ----------- //
		if (delay != 0)
		{
			RLoc = SampleCount - (int32_t) delay;
			append_to_file("output.txt",RLoc);
			printf("%d beat detected:\n", RLoc);
			//beat_mask[RLoc]=1;
			++Rcount;
		}
		else
		{
			RLoc = 0;
		}

		// -------- Toolbox comes with many helper functions for debugging, see PanTompkins.c for more details ---------- //
		s1 = PT_get_LPFilter_output();
		s2 = PT_get_HPFilter_output();
		s3 = PT_get_DRFilter_output();
		s4 = PT_get_SQRFilter_output();
		s5 = PT_get_MVFilter_output();

		ThI1 = PT_get_ThI1_output();
		SPKI = PT_get_SKPI_output();
		NPKI = PT_get_NPKI_output();
		ThF1 = PT_get_ThF1_output();

//		if (verbosity)
//			printf("%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", c, s1, s2, s3, s4, s5, RLoc, ThI1, SPKI, NPKI, ThF1);

//		fprintf_s(fptr_out, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", c, s1, s2, s3, s4, s5, RLoc, ThI1, SPKI, NPKI, ThF1);

	}
	printf("%d beats detected\n", Rcount);

//    //BWHighPass* filter = create_bw_high_pass_filter(5, 360, 5);
//    BWBandPass* filter = create_bw_band_pass_filter(5, 360, 0.5,30);
//    BWBandStop* filter2 = create_bw_band_stop_filter(5, 360, 53,57);
//
//   for(int i = 0; i < SIGNAL_SIZE; i++){
//      clean_signal  [i]= bw_band_pass(filter, ecg_signal[i]);
//      //clean_signal  [i]= apply_filter(b,a,ecg_signal[i]);
//    }
//    free_bw_band_pass(filter);
//
//    for(int i = 0; i < SIGNAL_SIZE; i++){
//      clean_signal  [i]= bw_band_stop(filter2, clean_signal[i]);
//      //clean_signal  [i]= apply_filter(b,a,ecg_signal[i]);
//    }
//    free_bw_band_stop(filter2);





    // write_array_to_file("output_array.bin", beat_mask, read_count);

    // Process float_array as needed
    // Note: There are 1000 float values in float_array

    return 0;
}


int read_floats_from_file(const char *filename, float **array, int read_count) {
    FILE *file;

    // Open the file for reading
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }

    // Allocate memory for the array based on read_count
    *array = (float *)malloc(read_count * sizeof(float));
    if (*array == NULL) {
        perror("Error allocating memory");
        fclose(file);
        return -1;
    }

    // Read the specified number of float values from the file
    for (int i = 0; i < read_count; i++) {
        if (fscanf(file, "%f", &(*array)[i]) != 1) {
            perror("Error reading float value from file");
            fclose(file);
            return -1;
        }
    }

    // Close the file
    fclose(file);

    return 0;
}


