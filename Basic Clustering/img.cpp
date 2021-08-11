#include <stdint.h>
#include <bits/stdc++.h>

using namespace std;
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

int main() {
    int width, height, bpp;

    uint8_t* rgb_image = stbi_load("test.jpg", &width, &height, &bpp, 3);
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            cout << rgb_image[i*sizeof(uint8_t) + j] << " ";
        }
        cout << endl;
    }
    stbi_image_free(rgb_image);

    return 0;
}