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

// Function to calculate the value
// of y
// y = poly[0] + x*poly[1] + x^2*poly[2] + ...
int calculate_Y(int x, vector<int>& poly)
{
    // Initializing y
    int y = 0;
    int temp = 1;

    // Iterating through the array
    for (auto coeff : poly) {

        // Computing the value of y
        y = (y + (coeff * temp));
        temp = (temp * x);
    }
    return y;
}

// Function to perform the secret
// sharing algorithm and encode the
// given secret
void secret_sharing(int S, vector<pair<int, int> >& points,
                    int N, int K)
{
    // A vector to store the polynomial
    // coefficient of K-1 degree
    vector<int> poly(K);

    // Randomly choose K - 1 numbers but
    // not zero and poly[0] is the secret
    // create polynomial for this

    poly[0] = S;

    for (int j = 1; j < K; ++j) {
        int p = 0;
        while (p == 0)

            // To keep the random values
            // in range not too high
            // we are taking mod with a
            // prime number around 1000
            p = (rand() % 997);

        // This is to ensure we did not
        // create a polynomial consisting
        // of zeroes.
        poly[j] = p;
    }

    // Generating N points from the
    // polynomial we created
    for (int j = 1; j <= N; ++j) {
        int x = j;
        int y = calculate_Y(x, poly);

        // Points created on sharing
        points[j - 1] = { x, y };
    }
}

// This structure is used for fraction
// part handling multiplication
// and addition of fractiontion
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

// Function to generate the secret
// back from the given points
// This function will use Lagrange Basis Polynomial
// Instead of finding the complete Polynomial
// We only required the poly[0] as our secret code,
// thus we can get rid of x terms
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
                l = l * temp;
            }
        }
        ans = ans + l;
    }

    // Return the secret
    return ans.num;
}

// Function to encode and decode the
// given secret by using the above
// defined functions
void operation(int S, int N, int K)
{

    // Vector to store the points
    vector<pair<int, int> > points(N);

    // Sharing of secret Code in N parts
    secret_sharing(S, points, N, K);

    cout << "Secret is divided to " << N
         << " Parts - " << endl;

    for (int i = 0; i < N; ++i) {
        cout << points[i].first << " "
             << points[i].second << endl;
    }

    cout << "We can generate Secret from any of "
         << K << " Parts" << endl;

    // Input any M points from these
    // to get back our secret code.
    int M = K;

    // M can be greater than or equal to threshold but
    // for this example we are taking for threshold
    if (M < K) {
        cout << "Points are less than threshold "
             << K << " Points Required" << endl;
    }

    int* x = new int[M];
    int* y = new int[M];

    // Input M points you will get the secret
    // Let these points are first M points from
    // the N points which we shared above
    // We can take any M points

    for (int i = 0; i < M; ++i) {
        x[i] = points[i].first;
        y[i] = points[i].second;
    }

    // Get back our result again.
    cout << "Our Secret Code is : "
         << Generate_Secret(x, y, M) << endl;
}

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


void test_Scalar_threshold(void) {
//    embedded_pairing::bls12_381::PowersOfX alphax;
//    Scalar alpha;
//    random_zpstar(alphax, alpha, random_bytes); // 使用alphax函数生成一个随机的标量alpha
//    std::cout << "Value of alpha: " << static_cast<int>(alpha) << std::endl;
//    printf("%d \n", alpha);
//
//    operation(static_cast<int>(alpha), 4, 2);
}


void test_wkdibe_encrypt_decrypt(void) {
    MasterKey msk;
    setup(p, msk, 10, false, random_bytes, 0);
    // msk拆分成阈值(n, t)-threshold
//    int n = 3;
//    int t = 2;
//    split_master_key(msk, n, t, msks);
//    // 使用msk_i派生密钥sk1_i
    keygen(sk1, p, msk, attrs2, random_bytes);
    // 重构派生的密钥，即合并sk1_i为sk1
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

//    test_wkdibe_encrypt_decrypt_master<true>();
    test_wkdibe_encrypt_decrypt();
    test_Scalar_threshold();
//    test_wkdibe_qualifykey();
//    test_wkdibe_nondelegablekey();
//    test_wkdibe_adjust();
//    test_wkdibe_sign();
//    test_wkdibe_marshal<true>("Marshal Compressed");
//    test_wkdibe_marshal<false>("Marshal Uncompressed");
//    printf("DONE\n");
}


