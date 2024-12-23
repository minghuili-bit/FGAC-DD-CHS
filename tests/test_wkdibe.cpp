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

#include <chrono>
#include <iostream>
#include <vector>

using namespace embedded_pairing::wkdibe;
using namespace std;

std::vector<int> lagrange_coefficients;

//using embedded_pairing::core::BigInt;
extern "C" {
    void random_bytes(void* buffer, size_t len);
    uint64_t current_time_nanos(void);
}

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

SecretKey SecretKeyPartial[4];


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

void test_Scalar_threshold(void) {
//    embedded_pairing::bls12_381::PowersOfX alphax;
//    Scalar alpha;
//    random_zpstar(alphax, alpha, random_bytes); // 使用alphax函数生成一个随机的标量alpha
//    std::cout << "Value of alpha: " << static_cast<int>(alpha) << std::endl;
//    printf("%d \n", alpha);
//
//    operation(static_cast<int>(alpha), 4, 2);

    Scalar scalar;
    scalar.std_words[0] = -12345;
//    scalar.std_words[0] = static_cast<uint32_t>(i);
    std::cout << "Value of alpha: " << static_cast<int>(scalar) << std::endl;

}

void test_wkdibe_encrypt_decrypt(void) {
    MasterKey msk;
    setup(p, msk, 10, false, random_bytes);

    std::cout << "Back Value: " << std::endl;

    for (int i = 0; i < 4; ++i) {
        cout << msk.points[i].first << " "
             << msk.points[i].second << endl;
    }

    int* x = new int[2];
    int* y = new int[2];

    for (int i = 0; i < 2; ++i) {
        x[i] = msk.points[i].first;
        y[i] = msk.points[i].second;
    }

    int revert_value = Generate_Secret(x, y, 2);

    std::cout << "revert_value: " << revert_value << std::endl;

    for (int i = 0; i < 2; ++i) {
        cout << lagrange_coefficients[i] << endl;
    }

    G1 result = G1::zero; // 初始化为零点
    for (int i = 0; i < 2; ++i) {

        Scalar scalar;
        scalar.std_words[0] = lagrange_coefficients[i];
        G1 temp;
        temp.multiply(msk.g2alpha1[i], scalar);
        result.add(result, temp);
//        result.add(temp);
//        msk.g2alpha1[i].multiply(scalar);
        G1 new_result;
        new_result.add(result, temp); // 使用基类的 add 函数
        result = new_result; // 更新 result
    }
    Scalar scalar1;
    scalar1.std_words[0] = revert_value;
    msk.g2alpha.multiply(p.g2, scalar1);

    if (embedded_pairing::bls12_381::Projective<embedded_pairing::bls12_381::Fq>::equal(result, msk.g2alpha)) {
        std::cout << "The results are identical." << std::endl;
    } else {
        std::cout << "The results differ." << std::endl;
    }

    keygen1(sk1, p, result, attrs2, random_bytes);
    // 恢复密钥，从msk.g2alpha1[i]恢复出来主密钥
//    for (int i = 0; i < 4; ++i) {
//        keygen1(sk1, p, msk.g2alpha1[i], attrs2, random_bytes);
//        SecretKeyPartial[i] = sk1;
//    }

    // 重构派生的密钥，即合并SecretKeyPartial[i]为sk1，难道是使用和g2乘积之后的结果来恢复，那太大了吧，即使能恢复

    GT msg;
    msg.random(random_bytes);

    Ciphertext c;
    encrypt(c, msg, p, attrs2, random_bytes);

    GT decrypted;
    decrypt(decrypted, c, sk1);

    if (GT::equal(msg, decrypted)) {
        printf("Decrypt Master: PASS\n");
    } else {
        printf("Decrypt Master: FAIL (original/decrypted messages differ)\n");
    }
}

extern "C" {
    void run_wkdibe_tests(void);
}

void run_wkdibe_tests() {
    init_test_wkdibe();

//    test_Scalar_threshold();
//    test_wkdibe_encrypt_decrypt_master<true>();
    test_wkdibe_encrypt_decrypt();
//    test_wkdibe_qualifykey();
//    test_wkdibe_nondelegablekey();
//    test_wkdibe_adjust();
//    test_wkdibe_sign();
//    test_wkdibe_marshal<true>("Marshal Compressed");
//    test_wkdibe_marshal<false>("Marshal Uncompressed");
//    printf("DONE\n");
}


