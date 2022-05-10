#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    if (x < 0) {
        x = 0;
    }
    if (x > im.w-1) {
        x = im.w-1;
    }
    if (y < 0) {
        y = 0;
    }
    if (y > im.h-1) {
        y = im.h-1;
    }
    if (c < 0) {
        c = 0;
    }
    if (c > im.c-1) {
        c = im.c-1;
    }
    return im.data[im.w * im.h * c + im.w * y + x];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    if (x < 0 || x > im.w-1 || y < 0 || y > im.h-1 || c < 0 || c > im.c-1) {
        return;
    } else {
        return im.data[im.w * im.h * c + im.w * y + x] = v;
    }
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    memcpy(copy.data, im.data, im.w * im.h * im.c * sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float R = get_pixel(im, i, j, 0);
            float G = get_pixel(im, i, j, 1);
            float B = get_pixel(im, i, j, 2);
            float Y = 0.299 * R + 0.587 * G + 0.114 * B;
            set_pixel(gray, i, j, 0, Y);
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            set_pixel(im, i, j, c, get_pixel(im, i, j, c) + v);
        }
    }
}

void clamp_image(image im)
{
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            for (int c = 0; c < im.c; c++) {
                float v = get_pixel(im, i, j, c);
                if (v < 0) {
                    v = 0.0;
                }
                if (v > 1) {
                    v = 1.0;
                }
                set_pixel(im, i, j, c, v);
            }
        }
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float R = get_pixel(im, i, j, 0);
            float G = get_pixel(im, i, j, 1);
            float B = get_pixel(im, i, j, 2);
            float V = three_way_max(R, G, B);
            float m = three_way_min(R,G,B);
            float C = V - m;
            float S = (R == 0.0 && G == 0.0 & B == 0.0) ? 0.0 : (C / V);
            float H_prime = 0.0;
            if (V == R) {
                H_prime = (G - B) / C;
            }
            if (V == G) {
                H_prime = (B - R) / C + 2;
            }
            if (V == B) {
                H_prime = (R - G) / C + 4;
            }
            float H = (H_prime < 0) ? (H_prime / 6) + 1 : (H_prime / 6);
            if (C == 0.0) {
                H = 0.0;
            }
            set_pixel(im, i, j, 0, H);
            set_pixel(im, i, j, 1, S);
            set_pixel(im, i, j, 2, V);
        }
    }
}

void hsv_to_rgb(image im)
{
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++) {
            float H = get_pixel(im, i, j, 0);
            float S = get_pixel(im, i, j, 1);
            float V = get_pixel(im, i, j, 2);
            H = H * 6;
            float Hi = floor(H);
            float F = H - Hi;
            float P = V * (1 - S);
            float Q = V * (1 - F * S);
            float T = V * (1 - (1 - F) * S);
            if (Hi == 0) {
                set_pixel(im, i, j, 0, V);
                set_pixel(im, i, j, 1, T);
                set_pixel(im, i, j, 2, P);
            } else if (Hi == 1) {
                set_pixel(im, i, j, 0, Q);
                set_pixel(im, i, j, 1, V);
                set_pixel(im, i, j, 2, P);
            } else if (Hi == 2) {
                set_pixel(im, i, j, 0, P);
                set_pixel(im, i, j, 1, V);
                set_pixel(im, i, j, 2, T);
            } else if (Hi == 3) {
                set_pixel(im, i, j, 0, P);
                set_pixel(im, i, j, 1, Q);
                set_pixel(im, i, j, 2, V);
            } else if (Hi == 4) {
                set_pixel(im, i, j, 0, T);
                set_pixel(im, i, j, 1, P);
                set_pixel(im, i, j, 2, V);
            } else if (Hi == 5) {
                set_pixel(im, i, j, 0, V);
                set_pixel(im, i, j, 1, P);
                set_pixel(im, i, j, 2, Q);
            }
        }
    }
}
