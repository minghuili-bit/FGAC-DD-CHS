#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "../include/bls12_381/fr.hpp"
#include "../include/bls12_381/fq.hpp"
#include "../include/bls12_381/fq2.hpp"
#include "../include/bls12_381/fq6.hpp"
#include "../include/bls12_381/fq12.hpp"
#include "../include/bls12_381/curve.hpp"
#include "../include/bls12_381/pairing.hpp"
#include "../include/wkdibe/api.hpp"
#include "../include/core/bigint.hpp"

#include <chrono>
#include <iostream>
#include <vector>

using namespace embedded_pairing::wkdibe;
using namespace std;
using namespace embedded_pairing::core;

std::vector<int> lagrange_coefficients;

//using embedded_pairing::core::BigInt;
extern "C" {
    void random_bytes(void* buffer, size_t len);
    uint64_t current_time_nanos(void);
}
const int threshold = 7;
const int totalShares = 11;

// 定义分数结构
struct Fraction {
    int numerator;   // 分子
    int denominator; // 分母
    Fraction(int num, int denom) : numerator(num), denominator(denom) {}

};



G1 harr[10];
Params p;

Attribute attr1arr[] = {{{.std_words = {15}}, 5, false}};
AttributeList attrs1;

Attribute attr2arr[] = {{{.std_words = {12}}, 3, false}, {{.std_words = {15}}, 5, false}};
AttributeList attrs2;

Attribute attr3arr[] = {{{.std_words = {12}}, 3, false}, {{.std_words = {7}}, 4, false}, {{.std_words = {15}}, 5, false}};
AttributeList attrs3;

FreeSlot b1arr[10];
SecretKey sk1;

FreeSlot b2arr[10];
SecretKey sk2;

FreeSlot b3arr[10];
SecretKey sk3;



void init_test_wkdibe(void) {
    p.h = harr;

    attrs1.length = 1;
    attrs1.attrs = attr1arr;
    attrs1.omitAllFromKeysUnlessPresent = false;

    attrs2.length = 2;
    attrs2.attrs = attr2arr;
    attrs2.omitAllFromKeysUnlessPresent = false;

    attrs3.length = 3;
    attrs3.attrs = attr3arr;
    attrs2.omitAllFromKeysUnlessPresent = false;

    sk1.b = b1arr;
    sk2.b = b2arr;
    sk3.b = b3arr;
}

/* Allows us to pass brace-enclosed initializer lists to macros. */
#define ARR(...) __VA_ARGS__

struct fraction {
    int num, den;

    // A fraction consists of a
    // numerator and a denominator
    fraction(int n, int d)
    {
        num = n, den = d;
    }

    // If the fraction is not
    // in its reduced form
    // reduce it by dividing
    // them with their GCD
    int gcd(int a, int b) {
        return b == 0 ? a : gcd(b, a % b);
    }

    void reduce_fraction(fraction& f) {
        int gcd_val = gcd(f.num, f.den);
        f.num /= gcd_val;
        f.den /= gcd_val;
    }

    // Performing multiplication on the
    // fraction
    fraction operator*(fraction f)
    {
        fraction temp(num * f.num, den * f.den);
        reduce_fraction(temp);
        return temp;
    }

    // Performing addition on the
    // fraction
    fraction operator+(fraction f)
    {
        fraction temp(num * f.den + den * f.num,
                      den * f.den);

        reduce_fraction(temp);
        return temp;
    }
};

int Generate_Secret(int x[], int y[], int M)
{

    fraction ans(0, 1);

    // Loop to iterate through the given
    // points
    for (int i = 0; i < M; ++i) {

        // Initializing the fraction
        fraction l(y[i], 1);
        for (int j = 0; j < M; ++j) {

            // Computing the lagrange terms
            if (i != j) {
                fraction temp(-x[j], x[i] - x[j]);

                lagrange_coefficients.push_back(temp.num / temp.den);

                l = l * temp;
            }
        }
        ans = ans + l;
    }

    // Return the secret
    return ans.num;
}

// GCD function to calculate greatest common divisor
int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// Function to compute Lagrange coefficient for index i
int computeLagrangeCoefficient(int i, const std::vector<int>& indices, int K) {
    int numerator = 1;   // Numerator of the Lagrange coefficient
    int denominator = 1; // Denominator of the Lagrange coefficient

    // We only need to consider the first K points, not the whole indices array
    for (size_t j = 0; j < K; ++j) {
        if (i != j) {
            numerator *= indices[j];                  // Multiply by x_j for the numerator
            denominator *= (indices[j] - indices[i]); // Multiply by (x_j - x_i) for the denominator

            // Simplify the fraction
            int gcd_value = gcd(numerator, denominator);
            numerator /= gcd_value;
            denominator /= gcd_value;
        }
    }

    // If denominator is not 1, throw an exception (shouldn't happen theoretically)
    if (denominator != 1) {
        throw std::runtime_error("Error: Denominator should be 1 after simplification");
    }

    return numerator;
}

void test_Scalar_threshold(void) {
//    embedded_pairing::bls12_381::PowersOfX alphax;
//    Scalar alpha;
//    random_zpstar(alphax, alpha, random_bytes); // 使用alphax函数生成一个随机的标量alpha
//    std::cout << "Value of alpha: " << static_cast<int>(alpha) << std::endl;
//    printf("%d \n", alpha);
//
//    operation(static_cast<int>(alpha), 4, 2);

    Scalar scalar;
    scalar.std_words[0] = -1234567890;
//    scalar.std_words[0] = static_cast<uint32_t>(i);
    std::cout << "Value of alpha: " << static_cast<int>(scalar) << std::endl;
    // 实现预计算功能
    // Initialize indices with values from 1 to 10
//    std::vector<int> indices;
//    for (int i = 0; i < 10; ++i) {
//        indices.push_back(i + 1); // Add i + 1 to the list of indices
//    }

    // The threshold value (K) is the number of lambdas you want to compute
//    int K = 6; // For example, compute 6 coefficients
//    std::vector<int> lambdas(K);

    // Compute the Lagrange coefficients for the first K indices
//    for (size_t i = 0; i < K; ++i) {
//        lambdas[i] = computeLagrangeCoefficient(i, indices, K);
//        std::cout << "Value of lambda[" << i+1 << "]: " << lambdas[i] << std::endl;
//    }


}
//void test_wkdibe_encrypt_decrypt(void) {
//    constexpr int num_iterations = 20; // 多次运行的次数
//    MasterKey msk;
//
//    auto start2 = std::chrono::high_resolution_clock::now();
//    ThresholdSetup(p, msk, 10, false, random_bytes, totalShares, threshold);
//    // 获取结束时间点
//    auto end2 = std::chrono::high_resolution_clock::now();
//    // 计算耗时，单位可以是微秒、毫秒或秒
//    std::chrono::duration<double, std::milli> elapsed2 = end2 - start2;
//    // 输出运行时间
//    std::cout << "setup运行时间: " << elapsed2.count() << " 毫秒" << std::endl;
//
//    auto start1 = std::chrono::high_resolution_clock::now();
//    std::vector<SecretKey> derivedPartialKeys(threshold);
//    for (int i = 0; i < threshold; ++i) {
//
//        DelKeyDer(sk1, p, msk.g2AlphaShares[i], attrs2, random_bytes);
//        derivedPartialKeys[i] = sk1;
//    }
//    // 获取结束时间点
//    auto end1 = std::chrono::high_resolution_clock::now();
//    // 计算耗时，单位可以是微秒、毫秒或秒
//    std::chrono::duration<double, std::milli> elapsed1 = end1 - start1;
//    // 输出运行时间
//    std::cout << "DelKeyDer运行时间: " << elapsed1.count() << " 毫秒" << std::endl;
//
//
//    auto start = std::chrono::high_resolution_clock::now();
//
////    DelKeyDer(sk1, p, msk.g2AlphaShares[i], attrs2, random_bytes);
//    SecretKey reconstructed_key;
////    Reconstruct( ,p, derivedPartialKeys, attrs2, random_bytes)
//    reconstructKey(reconstructed_key, p, derivedPartialKeys, threshold, attrs2, random_bytes);
//
//    // 计算拉格朗日系数
////    std::vector<int> indices;
////    for (int i = 0; i < threshold; ++i) {
////        indices.push_back(i + 1); // Add i + 1 to the list of indices
////    }
////
////    std::vector<int> lambdas(threshold);
////    for (size_t i = 0; i < threshold; ++i) {
//////        计算第i个拉格朗日系数，并将结果存储在 lambdas 数组中
////        lambdas[i] = computeLagrangeCoefficient(i, indices, threshold);
////        std::cout << "L_" << (i+1) << "(0) = " << lambdas[i]  << std::endl;
////    }
//    // 和a0无关的内容
////    reconstructed_key.a1 = derivedPartialKeys[0].a1;
////    reconstructed_key.l = derivedPartialKeys[0].l; // 假设所有部分密钥的长度一致
////    reconstructed_key.signatures = derivedPartialKeys[0].signatures;
////    reconstructed_key.bsig = derivedPartialKeys[0].bsig;
////    reconstructed_key.b = derivedPartialKeys[0].b;
//
//    // 开始a0的恢复过程
////    reconstructed_key.a0.copy(G1::zero);
//
//
//
//    // 遍历每个部分密钥
////    for (int i = 0; i < threshold; ++i) {
////
////        // 计算 λ_i * partial_a0s[i]
////        G1 temp;
////        Scalar scalar;
////        scalar.std_words[0] = lambdas[i];
////        temp.multiply(derivedPartialKeys[i].a0, scalar); // temp = λ_i * partial_a0s[i]
////        // 累加到结果中
////        #pragma omp critical
////        {
////            reconstructed_key.a0.add(reconstructed_key.a0, temp);
////        }
////    }
//
//    // 获取结束时间点
//    auto end = std::chrono::high_resolution_clock::now();
//    // 计算耗时，单位可以是微秒、毫秒或秒
//    std::chrono::duration<double, std::milli> elapsed = end - start;
//    // 输出运行时间
//    std::cout << "Reconstruct运行时间: " << elapsed.count() << " 毫秒" << std::endl;
//
//
//    GT msg;
//    msg.random(random_bytes);
//
//    Ciphertext c;
//    encrypt(c, msg, p, attrs2, random_bytes);
//
//    GT decrypted;
//    decrypt(decrypted, c, reconstructed_key);
//
//    if (GT::equal(msg, decrypted)) {
//        printf("Decrypt Master: PASS\n");
//    } else {
//        printf("Decrypt Master: FAIL (original/decrypted messages differ)\n");
//    }
//}
void test_wkdibe_encrypt_decrypt(void) {
    MasterKey msk;
    threshold_setup(p, msk, 10, false, random_bytes, totalShares, threshold);

    std::vector<SecretKey> derivedPartialKeys(threshold);
    for (int i = 0; i < threshold; ++i) {

        delegate_key_derive(sk1, p, msk.g2AlphaShares[i], attrs2, random_bytes);
        derivedPartialKeys[i] = sk1;
    }

    SecretKey reconstructedKey;
    reconstruct_key(reconstructedKey, p, derivedPartialKeys, threshold, attrs2, random_bytes);

    GT msg;
    msg.random(random_bytes);

    Ciphertext c;
    encrypt(c, msg, p, attrs2, random_bytes);

    GT decrypted;
    decrypt(decrypted, c, reconstructedKey);

    if (GT::equal(msg, decrypted)) {
        printf("Decrypt Master: PASS\n");
    } else {
        printf("Decrypt Master: FAIL (original/decrypted messages differ)\n");
    }
}

void test_wkdibe_encrypt_decrypt1(void) {
    using namespace std::chrono; // 简化命名

    MasterKey msk;

    // 统计 setup 的运行时间
    auto start = high_resolution_clock::now();
    setup(p, msk, 10, false, random_bytes);
    auto end = high_resolution_clock::now();
    printf("setup: %f ms\n", duration<double, std::milli>(end - start).count());

    // 统计 keygen 的运行时间
    start = high_resolution_clock::now();
    keygen(sk1, p, msk, attrs1, random_bytes);
    end = high_resolution_clock::now();
    printf("keygen: %f ms\n", duration<double, std::milli>(end - start).count());

    // 统计 qualifykey 的运行时间
    start = high_resolution_clock::now();
    qualifykey(sk2, p, sk1, attrs2, random_bytes);
    end = high_resolution_clock::now();
    printf("qualifykey: %f ms\n", duration<double, std::milli>(end - start).count());

    // 统计 msg.random 的运行时间
    GT msg;
    start = high_resolution_clock::now();
    msg.random(random_bytes);
    end = high_resolution_clock::now();
    printf("msg.random: %f ms\n", duration<double, std::milli>(end - start).count());

    // 统计 encrypt 的运行时间
    Ciphertext c;
    start = high_resolution_clock::now();
    encrypt(c, msg, p, attrs2, random_bytes);
    end = high_resolution_clock::now();
    printf("encrypt: %f ms\n", duration<double, std::milli>(end - start).count());

    // 统计 decrypt 的运行时间
    GT decrypted;
    start = high_resolution_clock::now();
    decrypt(decrypted, c, sk2);
    end = high_resolution_clock::now();
    printf("decrypt: %f ms\n", duration<double, std::milli>(end - start).count());

    // 比较结果
    start = high_resolution_clock::now();
    bool is_equal = GT::equal(msg, decrypted);
    end = high_resolution_clock::now();
    printf("GT::equal: %f ms\n", duration<double, std::milli>(end - start).count());

    if (is_equal) {
        printf("QualifyKey: PASS\n");
    } else {
        printf("QualifyKey: FAIL (original/decrypted messages differ)\n");
    }
}

extern "C" {
    void run_wkdibe_tests(void);
}

void run_wkdibe_tests() {
    init_test_wkdibe();
    test_wkdibe_encrypt_decrypt();

}


