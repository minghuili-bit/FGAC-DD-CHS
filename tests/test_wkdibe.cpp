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

using namespace embedded_pairing::wkdibe;
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

void test_wkdibe_encrypt_decrypt_master(void) {
    MasterKey msk;
    setup(p, msk, 10, false, random_bytes);

    GT msg;
    msg.random(random_bytes);

    Ciphertext c;
    encrypt(c, msg, p, attrs2, random_bytes);

    GT decrypted;
    decrypt_master(decrypted, c, msk);

    if (GT::equal(msg, decrypted)) {
        printf("Encrypt/Decrypt: PASS\n");
    } else {
        printf("Encrypt/Decrypt: FAIL (original/decrypted messages differ)\n");
    }
}

void test_wkdibe_encrypt_decrypt(void) {
    MasterKey msk;
    setup(p, msk, 10, false, random_bytes);
    keygen(sk1, p, msk, attrs2, random_bytes); // 根据master密钥派生sk1密钥

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

void test_wkdibe_qualifykey(void) {
    MasterKey msk;
    setup(p, msk, 10, false, random_bytes);
    keygen(sk1, p, msk, attrs1, random_bytes);
    qualifykey(sk2, p, sk1, attrs2, random_bytes);

    GT msg;
    msg.random(random_bytes);

    Ciphertext c;
    encrypt(c, msg, p, attrs2, random_bytes);

    GT decrypted;
    decrypt(decrypted, c, sk2);

    if (GT::equal(msg, decrypted)) {
        printf("QualifyKey: PASS\n");
    } else {
        printf("QualifyKey: FAIL (original/decrypted messages differ)\n");
    }
}

void test_wkdibe_nondelegablekey(void) {
    MasterKey msk;
    setup(p, msk, 10, false, random_bytes);
    nondelegable_keygen(sk1, p, msk, attrs1);
    nondelegable_qualifykey(sk2, p, sk1, attrs2);
    qualifykey(sk3, p, sk2, attrs3, random_bytes);

    GT msg;
    msg.random(random_bytes);

    Ciphertext c;
    encrypt(c, msg, p, attrs3, random_bytes);

    GT decrypted;
    decrypt(decrypted, c, sk3);

    if (GT::equal(msg, decrypted)) {
        printf("NonDelegableKey: PASS\n");
    } else {
        printf("NonDelegableKey: FAIL (original/decrypted messages differ)\n");
    }
}

void test_wkdibe_adjust(void) {
    auto start = std::chrono::high_resolution_clock::now();
    MasterKey msk;
    setup(p, msk, 10, false, random_bytes);
    keygen(sk1, p, msk, attrs1, random_bytes);
    nondelegable_qualifykey(sk2, p, sk1, attrs2);
    adjust_nondelegable(sk2, sk1, attrs2, attrs3);

    Precomputed precomputed;
    precompute(precomputed, p, attrs1);
    adjust_precomputed(precomputed, p, attrs1, attrs2);
    adjust_precomputed(precomputed, p, attrs2, attrs3);

    GT msg;
    msg.random(random_bytes);

    Ciphertext c;
    encrypt_precomputed(c, msg, p, precomputed, random_bytes);

    GT decrypted;
    decrypt(decrypted, c, sk2);

    if (GT::equal(msg, decrypted)) {
        printf("Adjust: PASS\n");
    } else {
        printf("Adjust: FAIL (original/decrypted messages differ)\n");
    }
    // 获取结束时间点
    auto end = std::chrono::high_resolution_clock::now();

    // 计算耗时（以微秒为单位）
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "运行时间: " << duration.count() << " 微秒" << std::endl;
}

void test_wkdibe_sign(void) {
    MasterKey msk;
    setup(p, msk, 10, true, random_bytes);
    keygen(sk1, p, msk, attrs2, random_bytes);

    Scalar msg;
    random_zpstar(msg, random_bytes);

    Signature s;
    sign(s, p, sk1, &attrs3, msg, random_bytes); // 属性2的密钥sk1签名属性3对应的信息

    if (verify(p, attrs3, s, msg)) {  // 签名验证：属性3+信息，用来验证签名s
        printf("Sign/Verify: PASS\n");
    } else {
        printf("Sign/Verify: FAIL (valid signature marked invalid)\n");
    }
}

template <bool compressed>
void test_wkdibe_marshal(const char* name) {
    {
        MasterKey msk;
        setup(p, msk, 10, false, random_bytes);

        {
            size_t pbuflen = p.getMarshalledLength<compressed>();
            uint8_t pbuf[pbuflen];
            p.marshal<compressed>(pbuf);
            if (p.setLength<compressed>(pbuf, pbuflen) == -1) {
                printf("%s: FAIL (could not set length for params)\n", name);
                return;
            }
            if (!p.unmarshal<compressed>(pbuf, true)) {
                printf("%s: FAIL (could not unmarshal params)\n", name);
                return;
            }
        }

        {
            uint8_t mskbuf[MasterKey::marshalledLength<compressed>];
            msk.marshal<compressed>(mskbuf);
            if (!msk.unmarshal<compressed>(mskbuf, true)) {
                printf("%s: FAIL (could not unmarshal master secret key)\n", name);
                return;
            }
        }

        keygen(sk1, p, msk, attrs1, random_bytes);

        {
            size_t sk1buflen = sk1.getMarshalledLength<compressed>();
            uint8_t sk1buf[sk1buflen];
            sk1.marshal<compressed>(sk1buf);
            if (sk1.setLength<compressed>(sk1buf, sk1buflen) == -1) {
                printf("%s: FAIL (could not set length for sk1)\n", name);
                return;
            }
            if (!sk1.unmarshal<compressed>(sk1buf, true)) {
                printf("%s: FAIL (could not unmarshal sk1)\n", name);
                return;
            }
        }

        qualifykey(sk2, p, sk1, attrs2, random_bytes);

        {
            size_t sk2buflen = sk2.getMarshalledLength<compressed>();
            uint8_t sk2buf[sk2buflen];
            sk2.marshal<compressed>(sk2buf);
            if (sk2.setLength<compressed>(sk2buf, sk2buflen) == -1) {
                printf("%s: FAIL (could not set length for sk2)\n", name);
                return;
            }
            if (!sk2.unmarshal<compressed>(sk2buf, true)) {
                printf("%s: FAIL (could not unmarshal sk2)\n", name);
                return;
            }
        }
    }

    {
        GT msg;
        msg.random(random_bytes);

        Ciphertext c;
        encrypt(c, msg, p, attrs2, random_bytes);

        {
            uint8_t cbuf[Ciphertext::marshalledLength<compressed>];
            c.marshal<compressed>(cbuf);
            if (!c.unmarshal<compressed>(cbuf, true)) {
                printf("%s: FAIL (could not unmarshal ciphertext)\n", name);
                return;
            }
        }

        GT decrypted;
        decrypt(decrypted, c, sk2);

        if (!GT::equal(msg, decrypted)) {
            printf("%s: FAIL (original/decrypted messages differ)\n", name);
            return;
        }
    }

    {
        Scalar msg2;
        random_zpstar(msg2, random_bytes);

        Signature s;
        sign(s, p, sk1, &attrs3, msg2, random_bytes);

        {
            uint8_t sbuf[Signature::marshalledLength<compressed>];
            s.marshal<compressed>(sbuf);
            if (!s.unmarshal<compressed>(sbuf, true)) {
                printf("%s: FAIL (could not unmarshal ciphertext)\n", name);
                return;
            }
        }

        if (!verify(p, attrs3, s, msg2)) {
            printf("%s: FAIL (valid signature marked invalid)\n", name);
            return;
        }
    }

    printf("%s: PASS\n", name);
}

extern "C" {
    void run_wkdibe_tests(void);
}

void run_wkdibe_tests() {
    init_test_wkdibe();

    test_wkdibe_encrypt_decrypt_master();
//    test_wkdibe_encrypt_decrypt();
//    test_wkdibe_qualifykey();
//    test_wkdibe_nondelegablekey();
//    test_wkdibe_adjust();
//    test_wkdibe_sign();
//    test_wkdibe_marshal<true>("Marshal Compressed");
//    test_wkdibe_marshal<false>("Marshal Uncompressed");
//    printf("DONE\n");
}


