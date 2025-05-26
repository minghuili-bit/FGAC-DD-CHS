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

#include "../include/lqibe/api.hpp"
#include "../include/lqibe/lqibe.h"

#include <chrono>
#include <iostream>
#include <vector>

using namespace std;
using namespace embedded_pairing::core;
using namespace embedded_pairing::wkdibe;
using namespace embedded_pairing::bls12_381;


extern "C" {
    void random_bytes(void* buffer, size_t len);
    uint64_t current_time_nanos(void);
}

void test_keygen(void) {
    G1 harr[10];
    Params params;
    params.h = harr;
    PowersOfX alphax;
    Scalar alpha;
    random_zpstar(alphax, alpha, random_bytes);
    params.g.random_generator(random_bytes);
    for(int i=0; i<8; i++){
        params.g1.multiply_frobenius(params.g, alphax);
    }
}

void test_dataEnc() {
    const size_t block_count = 100;
    const size_t block_size = 1024;
    const size_t key_len = 32;

    uint8_t seed[16];
    random_bytes(seed, sizeof(seed));

    bool all_pass = true;
    auto start_total = std::chrono::high_resolution_clock::now();
    double decrypt_total_ms = 0.0;

    for (size_t i = 0; i < block_count; ++i) {
        uint8_t data[block_size];
        random_bytes(data, sizeof(data));

        uint8_t key[key_len];
        random_bytes(key, sizeof(key));

        uint8_t ciphertext[block_size + 32];
        uint8_t recovered[block_size + 32];

        size_t enc_len = SymmetricCipher::encrypt(ciphertext, data, sizeof(data), key, key_len);
        size_t dec_len = SymmetricCipher::decrypt(recovered, ciphertext, enc_len, key, key_len);
    }
}


void test_tokenEnc() {
    Token t;
    t.id = 1001;
    t.prg = "G";
    t.uid = "user123";
    t.height = 5;
    t.offset = 3;
    t.seed = std::vector<uint8_t>(16);
    t.auth = std::vector<uint8_t>(32);
    random_bytes(t.seed.data(), t.seed.size());
    random_bytes(t.auth.data(), t.auth.size());
    t.gen_time = current_time_nanos();
    t.exp_time = t.gen_time + 60000000;  // 1 分钟后过期


    G1 harr[10];
    Params p;
    p.h = harr;

    Attribute attr2arr[] = {{{.std_words = {12}}, 3, false}, {{.std_words = {15}}, 5, false}};
    AttributeList attrs2;
    attrs2.length = 2;
    attrs2.attrs = attr2arr;
    attrs2.omitAllFromKeysUnlessPresent = false;

//    GT msg;
//    msg.random(random_bytes);
    Ciphertext c;
    std::vector<uint8_t> ciphertext;

    TokenEncryptor::encrypt_token(c, ciphertext, t, p, attrs2, random_bytes);


    FreeSlot b1arr[10];
    SecretKey sk1;
    sk1.b = b1arr;

    MasterKey msk;
    setup(p, msk, 10, false, random_bytes);
    keygen(sk1, p, msk, attrs2, random_bytes);

//    Token recovered = TokenEncryptor::decrypt_token(c, ciphertext, sk1);


}

void test_PerDel() {
    int threshold = 3;
    int totalShares = 5;
    MasterKey msk;
    G1 harr[10];
    Params p;
    p.h = harr;
    threshold_setup(p, msk, 10, false, random_bytes, totalShares, threshold);
}


void test_ThAuth() {
    int threshold = 7;
    int totalShares = 11;
    MasterKey msk;
    G1 harr[10];
    Params p;
    p.h = harr;
    Attribute attr2arr[] = {{{.std_words = {12}}, 3, false}, {{.std_words = {15}}, 5, false}};
    AttributeList attrs2;
    attrs2.length = 2;
    attrs2.attrs = attr2arr;
    attrs2.omitAllFromKeysUnlessPresent = false;
    FreeSlot b1arr[10];
    SecretKey sk1;
    sk1.b = b1arr;
    threshold_setup(p, msk, 10, false, random_bytes, totalShares, threshold);

    std::vector<SecretKey> derivedPartialKeys(threshold);
    for (int i = 0; i < threshold; ++i) {

        delegate_key_derive(sk1, p, msk.g2AlphaShares[i], attrs2, random_bytes);
        derivedPartialKeys[i] = sk1;
    }
    SecretKey reconstructedKey;
    reconstruct_key(reconstructedKey, p, derivedPartialKeys, threshold, attrs2, random_bytes);
}

extern "C" {
    void run_integration_tests(void);
}

void run_integration_tests() {
    test_keygen();
    test_dataEnc();
    test_tokenEnc();
    test_PerDel();
    test_ThAuth();
}

