#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include "Scene.hpp"
#include "Renderer.hpp"
#include <immintrin.h>


inline float deg2rad(const float& deg) { return deg * M_PI / 180.0; }

const float EPSILON = 0.00001;

// The main render function. This where we iterate over all pixels in the image,
// generate primary rays and cast these rays into the scene. The content of the
// framebuffer is saved to a file.
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = tan(deg2rad(scene.fov * 0.5));
    float imageAspectRatio = scene.width / (float)scene.height;
    Vector3f eye_pos(278, 273, -800);
    int m = 0;

    // change the spp value to change sample ammount
    int spp = 4;
    std::cout << "SPP: " << spp << "\n";
    for (uint32_t j = 0; j < scene.height; ++j) {
        for (uint32_t i = 0; i < scene.width; ++i) {
            // generate primary ray direction
            float x = (2 * (i + 0.5) / (float)scene.width - 1) *
                      imageAspectRatio * scale;
            float y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

            Vector3f dir = normalize(Vector3f(-x, y, 1));
            for (int k = 0; k < spp; k++){
                framebuffer[m] += scene.castRay(Ray(eye_pos, dir), 0) / spp;  
            }
            m++;
        }
        UpdateProgress(j / (float)scene.height);
    }
    //for (uint32_t j = 0; j < scene.height; j += 8) {
    //    for (uint32_t i = 0; i < scene.width; i += 8) {
    //        // generate primary ray direction
    //        __m256i idx = _mm256_set_epi32(i + 7, i + 6, i + 5, i + 4, i + 3, i + 2, i + 1, i);
    //        __m256i jdx = _mm256_set_epi32(j + 7, j + 6, j + 5, j + 4, j + 3, j + 2, j + 1, j);
    //        __m256i widthdx = _mm256_set1_epi32(scene.width);
    //        __m256i heightdx = _mm256_set1_epi32(scene.height);
    //        __m256 scalex = _mm256_set1_ps(scale);
    //        __m256 imageAspectRatiox = _mm256_set1_ps(imageAspectRatio);
    //        __m256 halfx = _mm256_set1_ps(0.5f);
    //        __m256 onex = _mm256_set1_ps(1.0f);
    //        __m256 monex = _mm256_set1_ps(-1.0f);
    //        __m256 twox = _mm256_set1_ps(2.0f);
    //        __m256 fx = _mm256_cvtepi32_ps(idx);
    //        __m256 fy = _mm256_cvtepi32_ps(jdx);
    //        fx = _mm256_mul_ps(_mm256_add_ps(fx, halfx), twox);//2*(i+0.5)
    //        fy = _mm256_mul_ps(_mm256_add_ps(fy, halfx), twox);//2*(j+0.5)
    //        fx = _mm256_div_ps(fx, _mm256_cvtepi32_ps(widthdx));//2*(i+0.5)/(float)scene.width
    //        fy = _mm256_div_ps(fy, _mm256_cvtepi32_ps(heightdx));//2 * (j + 0.5) / (float)scene.height
    //        fx = _mm256_sub_ps(fx, onex); //2*(i+0.5)/(float)scene.width-1
    //        fy = _mm256_add_ps(_mm256_mul_ps(fy, monex), onex);//(1 - 2 * (j + 0.5) / (float)scene.height)
    //        fx = _mm256_mul_ps(fx, imageAspectRatiox);
    //        fx = _mm256_mul_ps(fx, scalex);
    //        fy = _mm256_mul_ps(fy, scalex);
    //        float vecx[8];
    //        float vecy[8];
    //        _mm256_store_ps(vecx, fx);
    //        _mm256_store_ps(vecy, fy);
    //        Vector3f dir[64];
    //        int q = 0;
    //        for (int t = 0; t < 8; t++) {
    //            for (int p = 0; p < 8; p++) {
    //                dir[q] = normalize(Vector3f(-vecx[p], vecy[t], 1));
    //                q++;
    //            }
    //        }
    //        int r[8] = { m,m + 1,m + 2,m + 3,m + 4,m + 5,m + 6,m + 7 };
    //        for (int k = 0; k < spp; k++) {
    //            framebuffer[r[0]]+= scene.castRay(Ray(eye_pos, dir[0]), 0) / spp;
    //            framebuffer[r[1]] += scene.castRay(Ray(eye_pos, dir[1]), 0) / spp;
    //            framebuffer[r[2]] += scene.castRay(Ray(eye_pos, dir[2]), 0) / spp;
    //            framebuffer[r[3]] += scene.castRay(Ray(eye_pos, dir[3]), 0) / spp;
    //            framebuffer[r[4]] += scene.castRay(Ray(eye_pos, dir[4]), 0) / spp;
    //            framebuffer[r[5]] += scene.castRay(Ray(eye_pos, dir[5]), 0) / spp;
    //            framebuffer[r[6]] += scene.castRay(Ray(eye_pos, dir[6]), 0) / spp;
    //            framebuffer[r[7]] += scene.castRay(Ray(eye_pos, dir[7]), 0) / spp;
    //            framebuffer[r[0]+scene.width] += scene.castRay(Ray(eye_pos, dir[8]), 0) / spp;
    //            framebuffer[r[1]+scene.width] += scene.castRay(Ray(eye_pos, dir[9]), 0) / spp;
    //            framebuffer[r[2]+scene.width] += scene.castRay(Ray(eye_pos, dir[10]), 0) / spp;
    //            framebuffer[r[3]+scene.width] += scene.castRay(Ray(eye_pos, dir[11]), 0) / spp;
    //            framebuffer[r[4]+scene.width] += scene.castRay(Ray(eye_pos, dir[12]), 0) / spp;
    //            framebuffer[r[5]+scene.width] += scene.castRay(Ray(eye_pos, dir[13]), 0) / spp;
    //            framebuffer[r[6]+scene.width] += scene.castRay(Ray(eye_pos, dir[14]), 0) / spp;
    //            framebuffer[r[7]+scene.width] += scene.castRay(Ray(eye_pos, dir[15]), 0) / spp;
    //            framebuffer[r[0] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[16]), 0) / spp;
    //            framebuffer[r[1] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[17]), 0) / spp;
    //            framebuffer[r[2] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[18]), 0) / spp;
    //            framebuffer[r[3] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[19]), 0) / spp;
    //            framebuffer[r[4] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[20]), 0) / spp;
    //            framebuffer[r[5] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[21]), 0) / spp;
    //            framebuffer[r[6] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[22]), 0) / spp;
    //            framebuffer[r[7] + 2*scene.width] += scene.castRay(Ray(eye_pos, dir[23]), 0) / spp;
    //            framebuffer[r[0] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[24]), 0) / spp;
    //            framebuffer[r[1] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[25]), 0) / spp;
    //            framebuffer[r[2] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[26]), 0) / spp;
    //            framebuffer[r[3] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[27]), 0) / spp;
    //            framebuffer[r[4] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[28]), 0) / spp;
    //            framebuffer[r[5] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[29]), 0) / spp;
    //            framebuffer[r[6] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[30]), 0) / spp;
    //            framebuffer[r[7] + 3 * scene.width] += scene.castRay(Ray(eye_pos, dir[31]), 0) / spp;
    //            framebuffer[r[0] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[32]), 0) / spp;
    //            framebuffer[r[1] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[33]), 0) / spp;
    //            framebuffer[r[2] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[34]), 0) / spp;
    //            framebuffer[r[3] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[35]), 0) / spp;
    //            framebuffer[r[4] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[36]), 0) / spp;
    //            framebuffer[r[5] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[37]), 0) / spp;
    //            framebuffer[r[6] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[38]), 0) / spp;
    //            framebuffer[r[7] + 4 * scene.width] += scene.castRay(Ray(eye_pos, dir[39]), 0) / spp;
    //            framebuffer[r[0] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[40]), 0) / spp;
    //            framebuffer[r[1] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[41]), 0) / spp;
    //            framebuffer[r[2] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[42]), 0) / spp;
    //            framebuffer[r[3] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[43]), 0) / spp;
    //            framebuffer[r[4] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[44]), 0) / spp;
    //            framebuffer[r[5] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[45]), 0) / spp;
    //            framebuffer[r[6] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[46]), 0) / spp;
    //            framebuffer[r[7] + 5 * scene.width] += scene.castRay(Ray(eye_pos, dir[47]), 0) / spp;
    //            framebuffer[r[0] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[48]), 0) / spp;
    //            framebuffer[r[1] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[49]), 0) / spp;
    //            framebuffer[r[2] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[50]), 0) / spp;
    //            framebuffer[r[3] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[51]), 0) / spp;
    //            framebuffer[r[4] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[52]), 0) / spp;
    //            framebuffer[r[5] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[53]), 0) / spp;
    //            framebuffer[r[6] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[54]), 0) / spp;
    //            framebuffer[r[7] + 6 * scene.width] += scene.castRay(Ray(eye_pos, dir[55]), 0) / spp;
    //            framebuffer[r[0] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[56]), 0) / spp;
    //            framebuffer[r[1] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[57]), 0) / spp;
    //            framebuffer[r[2] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[58]), 0) / spp;
    //            framebuffer[r[3] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[59]), 0) / spp;
    //            framebuffer[r[4] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[60]), 0) / spp;
    //            framebuffer[r[5] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[61]), 0) / spp;
    //            framebuffer[r[6] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[62]), 0) / spp;
    //            framebuffer[r[7] + 7 * scene.width] += scene.castRay(Ray(eye_pos, dir[63]), 0) / spp;
    //        }
    //        m += 8;
    //    }
    //    UpdateProgress(j / (float)scene.height);
    //    m += 7 * scene.width;
    //}
    UpdateProgress(1.f);

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].x), 0.6f));
        color[1] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].y), 0.6f));
        color[2] = (unsigned char)(255 * std::pow(clamp(0, 1, framebuffer[i].z), 0.6f));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
