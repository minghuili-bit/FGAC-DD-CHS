/*
 * Copyright (c) 2018, Sam Kumar <samkumar@cs.berkeley.edu>
 * Copyright (c) 2018, University of California, Berkeley
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef EMBEDDED_PAIRING_WKDIBE_API_HPP_
#define EMBEDDED_PAIRING_WKDIBE_API_HPP_

#include <stddef.h>
#include <stdint.h>
#include <vector>

#include "core/bigint.hpp"
#include "bls12_381/fr.hpp"
#include "bls12_381/fq12.hpp"
#include "bls12_381/curve.hpp"
#include "bls12_381/pairing.hpp"
#include "bls12_381/decomposition.hpp"

namespace embedded_pairing::wkdibe {


    typedef bls12_381::G1 G1;
    typedef bls12_381::G2 G2;
    typedef bls12_381::G1Affine G1Affine;
    typedef bls12_381::G2Affine G2Affine;
    typedef bls12_381::Fq12 GT;
    typedef core::BigInt<256> ID;
    typedef core::BigInt<256> Scalar;

    const core::BigInt<256> group_order = bls12_381::Fr::p_value;


    struct Attribute {
        ID id;
        uint32_t idx;
        bool omitFromKeys;
    };

    struct AttributeList {
        Attribute* attrs;
        size_t length;
        bool omitAllFromKeysUnlessPresent;
    };

    struct Params {
        G2 g;
        G2 g1;
        G1 g2;
        G1 g3;
        GT pairing;

        G1 hsig;
        bool signatures;

        G1* h;
        int l;

        template <bool compressed>
        void marshal(void* buffer) const;

        template <bool compressed>
        bool unmarshal(const void* buffer, bool checked);

        template <bool compressed>
        inline int setLength(const void* marshalled, size_t marshalledLength) {
            int len = Params::unmarshalledLength<compressed>(marshalled, marshalledLength);
            if (len != -1) {
                this->l = len;
            }
            return len;
        }

        template <bool compressed>
        inline size_t getMarshalledLength(void) const {
            return Params::marshalledLength<compressed>(this->l, this->signatures);
        }

        template <bool compressed>
        static constexpr size_t marshalledLengthMinimum = 1 + 2 * bls12_381::Encoding<G1Affine, compressed>::size + 2 * bls12_381::Encoding<G2Affine, compressed>::size;

        template <bool compressed>
        static constexpr int unmarshalledLength(const void* marshalled, size_t marshalledLength) {
            uint8_t firstByte = *((uint8_t*) marshalled);
            size_t withoutLength = Params::marshalledLengthMinimum<compressed> + (compressed ? 0 : sizeof(GT)) + (firstByte == 0 ? 0 : bls12_381::Encoding<G1Affine, compressed>::size);
            if (marshalledLength < withoutLength) {
                return -1;
            }
            size_t hsize = marshalledLength - withoutLength;
            return (hsize % bls12_381::Encoding<G1Affine, compressed>::size) == 0 ? (hsize / bls12_381::Encoding<G1Affine, compressed>::size) : -1;
        }

        template <bool compressed>
        static constexpr size_t marshalledLength(int length, bool signatures) {
            return Params::marshalledLengthMinimum<compressed> + (compressed ? 0 : sizeof(GT)) + ((signatures ? 1 : 0) + length) * bls12_381::Encoding<G1Affine, compressed>::size;
        }
    };

    struct Ciphertext {
        GT a;
        G2 b;
        G1 c;

        template <bool compressed>
        void marshal(void* buffer) const;

        template <bool compressed>
        bool unmarshal(const void* buffer, bool checked);

        template <bool compressed>
        static constexpr size_t marshalledLength = bls12_381::Encoding<G1Affine, compressed>::size + bls12_381::Encoding<G2Affine, compressed>::size + sizeof(GT);
    };

    struct Signature {
        G1 a0;
        G2 a1;

        template <bool compressed>
        void marshal(void* buffer) const;

        template <bool compressed>
        bool unmarshal(const void* buffer, bool checked);

        template <bool compressed>
        static constexpr size_t marshalledLength = bls12_381::Encoding<G1Affine, compressed>::size + bls12_381::Encoding<G2Affine, compressed>::size;
    };

    struct FreeSlot {
        G1 hexp;
        uint32_t idx;

        template <bool compressed>
        static constexpr size_t marshalledLength = 4 + bls12_381::Encoding<G1Affine, compressed>::size;

        template <bool compressed>
        void marshal(void* buffer) const;

        template <bool compressed>
        bool unmarshal(const void* buffer, bool checked);
    };

    struct SecretKey {
        G1 a0;
        G2 a1;

        int l;
        bool signatures;
        G1 bsig;
        FreeSlot* b;

        template <bool compressed>
        void marshal(void* buffer) const;

        template <bool compressed>
        bool unmarshal(const void* buffer, bool checked);

        template <bool compressed>
        int setLength(const void* marshalled, size_t marshalledLength) {
            int len = SecretKey::unmarshalledLength<compressed>(marshalled, marshalledLength);
            if (len != -1) {
                this->l = len;
            }
            return len;
        }

        template <bool compressed>
        inline size_t getMarshalledLength(void) const {
            return SecretKey::marshalledLength<compressed>(this->l, this->signatures);
        }

        template <bool compressed>
        static constexpr size_t marshalledLengthMinimum = 1 + bls12_381::Encoding<G1Affine, compressed>::size + bls12_381::Encoding<G2Affine, compressed>::size;

        template <bool compressed>
        static constexpr int unmarshalledLength(const void* marshalled, size_t marshalledLength) {
            uint8_t firstByte = *((uint8_t*) marshalled);
            size_t withoutLength = SecretKey::marshalledLengthMinimum<compressed> + (firstByte == 0 ? 0 : bls12_381::Encoding<G1Affine, compressed>::size);
            if (marshalledLength < withoutLength) {
                return -1;
            }
            size_t bsize = marshalledLength - withoutLength;
            return (bsize % FreeSlot::marshalledLength<compressed>) == 0 ? (bsize / FreeSlot::marshalledLength<compressed>) : -1;
        }

        template <bool compressed>
        static constexpr size_t marshalledLength(int length, bool signatures) {
            return SecretKey::marshalledLengthMinimum<compressed> + length * FreeSlot::marshalledLength<compressed> + (signatures ? 1 : 0) * bls12_381::Encoding<G1Affine, compressed>::size;
        }
    };

    struct MasterKey {
        G1 g2alpha;

        std::vector<G1> g2AlphaShares; // 动态数组

//        std::vector<std::pair<int, int> > points;

        template <bool compressed>
        void marshal(void* buffer) const;

        template <bool compressed>
        bool unmarshal(const void* buffer, bool checked);

        template <bool compressed>
        static constexpr size_t marshalledLength = bls12_381::Encoding<G1Affine, compressed>::size;
    };

    struct Precomputed {
        G1 prodexp;
    };

    struct Token {
        uint32_t id;
        std::string prg;
        std::string uid;
        uint32_t height;
        uint32_t offset;
        std::vector<uint8_t> seed;
        uint64_t gen_time;
        uint64_t exp_time;
        std::vector<uint8_t> auth;
    };

    inline void scalar_hash_reduce(Scalar& x) {
        bls12_381::Fr* target = reinterpret_cast<bls12_381::Fr*>(&x);
        target->hash_reduce();
    }

    inline void random_zpstar(Scalar& s, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::Fr* target = reinterpret_cast<bls12_381::Fr*>(&s);
        target->random(get_random_bytes);
    }

    inline void random_zpstar(bls12_381::PowersOfX& __restrict powers, Scalar& __restrict s, void (*get_random_bytes)(void*, size_t)) {
        powers.random(s, get_random_bytes);
    }

    inline void random_g1(G1& g1, void (*get_random_bytes)(void*, size_t)) {
        g1.random_generator(get_random_bytes);
    }

    inline void random_g2(G2& g2, void (*get_random_bytes)(void*, size_t)) {
        g2.random_generator(get_random_bytes);
    }

    inline void random_gt(GT& gt, void (*get_random_bytes)(void*, size_t)) {
        Scalar s;
        gt.random_gt(s, bls12_381::generator_pairing, get_random_bytes);
    }

    class KeyDerivationTree {
    public:
        static std::vector<std::vector<uint8_t>> generate_keys(const uint8_t* seed, size_t seed_len, size_t num_keys, size_t key_len);
    };

    class SymmetricCipher {
    public:
        static size_t encrypt(uint8_t* out, const uint8_t* data, size_t len, const uint8_t* key, size_t key_len);
        static size_t decrypt(uint8_t* out, const uint8_t* ciphertext, size_t len, const uint8_t* key, size_t key_len);
    };

    class TokenEncryptor {
    public:
        static std::vector <uint8_t> serialize(const Token &token);

        static Token deserialize(const std::vector <uint8_t> &buffer);

        static void encrypt_token(Ciphertext &out_ct,
                                  std::vector <uint8_t> &out_sym_cipher,
                                  const Token &token,
                                  const Params &p,
                                  const AttributeList &attrs,
                                  void (*random_bytes)(void *, size_t));

        static Token decrypt_token(const Ciphertext& ct,
                                   const std::vector<uint8_t>& sym_cipher,
                                   const SecretKey& sk);
    };
    void setup(Params &params, MasterKey &msk, int l, bool signatures, void (*get_random_bytes)(void *, size_t));
    void threshold_setup(Params &params, MasterKey &msk, int l, bool signatures, void (*get_random_bytes)(void *, size_t), int N, int M);
    void delegate_key_derive(SecretKey& sk, const Params& params, const G1& msk, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t));
    void reconstruct_key(SecretKey& reconstructed_key, const Params& params, const std::vector<SecretKey>& derivedPartialKeys, int threshold, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t));
    void keygen(SecretKey& sk, const Params& params, const MasterKey& msk, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t));
    void qualifykey(SecretKey& qualified, const Params& params, const SecretKey& sk, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t));
    void nondelegable_keygen(SecretKey& sk, const Params& params, const MasterKey& msk, const AttributeList& attrs);
    void nondelegable_qualifykey(SecretKey& qualified, const Params& params, const SecretKey& sk, const AttributeList& attrs);
    void adjust_nondelegable(SecretKey& sk, const SecretKey& parent, const AttributeList& from, const AttributeList& to);

    void precompute(Precomputed& precomputed, const Params& params, const AttributeList& attrs);
    void adjust_precomputed(Precomputed& precomputed, const Params& params, const AttributeList& from, const AttributeList& to);
    void resamplekey(SecretKey& resampled, const Params& params, const Precomputed& precomputed, const SecretKey& sk, bool supportFurtherQualification, void (*get_random_bytes)(void*, size_t));

    void encrypt(Ciphertext& ciphertext, const GT& message, const Params& params, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t));
    void encrypt_precomputed(Ciphertext& ciphertext, const GT& message, const Params& params, const Precomputed& precomputed, void (*get_random_bytes)(void*, size_t));
    void decrypt(GT& message, const Ciphertext& ciphertext, const SecretKey& sk);
    void decrypt_master(GT& message, const Ciphertext& ciphertext, const MasterKey& msk);

    void sign(Signature& signature, const Params& params, const SecretKey& sk, const AttributeList* attrs, const Scalar& message, void (*get_random_bytes)(void*, size_t));
    void sign_precomputed(Signature& signature, const Params& params, const SecretKey& sk, const AttributeList* attrs, const Precomputed& precomputed, const Scalar& message, void (*get_random_bytes)(void*, size_t));
    bool verify(const Params& params, const AttributeList& attrs, const Signature& signature, const Scalar& message);
    bool verify_precomputed(const Params& params, const Precomputed& precomputed, const Signature& signature, const Scalar& message);
//    std::vector<std::vector<uint8_t>> KeyDerivationTree::generate_keys(const uint8_t* seed, size_t seed_len, size_t num_keys, size_t key_len)
//    void SymmetricCipher::encrypt(uint8_t* out, const uint8_t* data, size_t len, const uint8_t* key, size_t key_len)
//    void SymmetricCipher::decrypt(uint8_t* out, const uint8_t* ciphertext, size_t len, const uint8_t* key, size_t key_len)



    }

#endif
