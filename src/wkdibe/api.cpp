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

#include "wkdibe/api.hpp"

#include <stddef.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include "bls12_381/pairing.hpp"
#include "bls12_381/wnaf.hpp"
#include "bls12_381/decomposition.hpp"
#include "../../include/wkdibe/api.hpp"
#include <openssl/sha.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <cstring>
#include <stdexcept>

using namespace std;

namespace embedded_pairing::wkdibe {

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

    void secret_sharing(int S, int N, int K, std::vector<Scalar>& shares)
{

        vector<int> poly(K);
        poly[0] = S;

        for (int j = 1; j < K; ++j) {
            int p = 0;
            while (p == 0)
                p = (rand() % 997);
            poly[j] = p;
        }

        // Generate N y-values of the polynomial at x = 1, 2, ..., N
        for (int j = 1; j <= N; ++j) {
            int y = calculate_Y(j, poly);
            Scalar s;
            s.std_words[0] = static_cast<uint32_t>(y);
            for (size_t i = 1; i < sizeof(s.std_words)/sizeof(uint32_t); ++i)
            s.std_words[i] = 0;
            shares[j - 1] = s;
}
    }


    void threshold_setup(Params& params, MasterKey& msk, int l, bool signatures, void (*get_random_bytes)(void*, size_t), int totalShares, int threshold) {
        bls12_381::PowersOfX alphax;
        // Sample random scalar alpha in Z_p
        Scalar alpha;
        random_zpstar(alphax, alpha, get_random_bytes);
        params.g.random_generator(get_random_bytes);
        params.g1.multiply_frobenius(params.g, alphax);
        // Initialize group generator g2
        params.g2.random_generator(get_random_bytes);

//        std::cout << "Value of alpha: " << static_cast<int>(alpha) << std::endl;

        // Split alpha via Shamir sharing into shares
        std::vector<Scalar> shares(totalShares);
        secret_sharing(static_cast<int>(alpha), totalShares, threshold, shares);

//        cout << "Secret is divided to " << totalShares
//             << " Parts - " << endl;

        msk.g2AlphaShares.resize(totalShares);
        // Compute and store g_2^shares[i]
        for (int i = 0; i < totalShares; i++) {
//            auto mark = Scalar{.std_words = {(uint32_t)shares[i].first}};
//            auto secret = Scalar{.std_words = {(uint32_t)shares[i]}};
//            std::cout << "Value of i: " << static_cast<int>(mark) << std::endl;
//            std::cout << "Value of alpha: " << static_cast<int>(shares[i]) << std::endl;
            msk.g2AlphaShares[i].multiply(params.g2, shares[i]);
        }


        params.g3.random_generator(get_random_bytes);

        G1Affine g2affine;
        G2Affine g1affine;
        g2affine.from_projective(params.g2);
        g1affine.from_projective(params.g1);
        bls12_381::pairing(params.pairing, g2affine, g1affine);

        params.l = l;
        params.signatures = signatures;
        if (signatures) {
            params.hsig.random_generator(get_random_bytes);
        } else {
            params.hsig.copy(G1::zero);
        }
        for (int i = 0; i != l; i++) {
            params.h[i].random_generator(get_random_bytes);
        }
    }

    void calculatePatternToKey(SecretKey& sk, const Params& params, const AttributeList& attrs,
                              Scalar r, int& out_j) {
        G1 temp;
        int j = 0, k = 0;

        for (int i = 0; i != params.l; i++) {
            if (k != attrs.length && attrs.attrs[k].idx == i) {
                if (!attrs.attrs[k].omitFromKeys) {
                    temp.multiply(params.h[i], attrs.attrs[k].id);
                    sk.a0.add(sk.a0, temp);
                }
                k++;
            } else if (!attrs.omitAllFromKeysUnlessPresent) {
                sk.b[j].idx = i;
                sk.b[j].hexp.multiply(params.h[i], r);
                j++;
            }
        }

        sk.l = j;
        out_j = j;  // Optional: return j if needed by caller
    }

    int gcd(int a, int b) {
        while (b != 0) {
            int temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }

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
    // Function to compute Lagrange coefficients for threshold shares
    std::vector<Scalar> computeLagrangeCoefficients(int threshold) {
        std::vector<int> indices;
        for (int i = 0; i < threshold; ++i) {
            indices.push_back(i + 1);  // Add index i+1 to the list
        }

        std::vector<Scalar> lambdas(threshold);
        for (size_t i = 0; i < threshold; ++i) {
            int coeff = computeLagrangeCoefficient(i, indices, threshold);  // Compute Lagrange coefficient for each index
            Scalar scalar;
            scalar.std_words[0] = coeff;
            lambdas[i] = scalar;
//            std::cout << "L_" << (i + 1) << "(0) = " << lambdas[i] << std::endl;  // Optional debug output
        }

        return lambdas;
    }

    void delegate_key_derive(SecretKey& derivedPartialKey, const Params& params, const G1& g2AlphaShares, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX rx;
        Scalar r;
        G1 temp;
        random_zpstar(rx, r, get_random_bytes);
        // 打印

        derivedPartialKey.a0.copy(params.g3);
        int j = 0;
        calculatePatternToKey(derivedPartialKey, params, attrs, r, j);
        derivedPartialKey.l = j;
        derivedPartialKey.signatures = params.signatures;
        if (derivedPartialKey.signatures) {
            derivedPartialKey.bsig.multiply(params.hsig, r);
        } else {
            derivedPartialKey.bsig.copy(G1::zero);
        }
        derivedPartialKey.a0.multiply(derivedPartialKey.a0, r);

        derivedPartialKey.a0.add(derivedPartialKey.a0, g2AlphaShares);
        derivedPartialKey.a1.multiply_frobenius(params.g, rx);
    }

    void reconstruct_key(SecretKey& reconstructedKey, const Params& params, const std::vector<SecretKey>& derivedPartialKeys, int threshold, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t)) {

        reconstructedKey.a0.copy(G1::zero);

        // Calculate the Lagrange coefficients for interpolation
        std::vector<Scalar> lambdas = computeLagrangeCoefficients(threshold);

        // Use Lagrange interpolation to recover a0
        for (int i = 0; i < threshold; ++i) {

            G1 temp;
//            Scalar scalar;
//            scalar.std_words[0] = lambdas[i];
            temp.multiply(derivedPartialKeys[i].a0, lambdas[i]); // temp = λ_i * partial_a0s[i]
            reconstructedKey.a0.add(reconstructedKey.a0, temp);
        }

        reconstructedKey.a1 = derivedPartialKeys[0].a1;
        reconstructedKey.l = derivedPartialKeys[0].l;
        reconstructedKey.signatures = derivedPartialKeys[0].signatures;
        reconstructedKey.bsig = derivedPartialKeys[0].bsig;
        reconstructedKey.b = derivedPartialKeys[0].b;


    }


    void setup(Params& params, MasterKey& msk, int l, bool signatures, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX alphax;
        Scalar alpha;
        random_zpstar(alphax, alpha, get_random_bytes);
        params.g.random_generator(get_random_bytes);
        params.g1.multiply_frobenius(params.g, alphax);
        params.g2.random_generator(get_random_bytes);
        msk.g2alpha.multiply(params.g2, alpha);
        params.g3.random_generator(get_random_bytes);

        G1Affine g2affine;
        G2Affine g1affine;
        g2affine.from_projective(params.g2);
        g1affine.from_projective(params.g1);
        bls12_381::pairing(params.pairing, g2affine, g1affine);

        params.l = l;
        params.signatures = signatures;
        if (signatures) {
            params.hsig.random_generator(get_random_bytes);
        } else {
            params.hsig.copy(G1::zero);
        }
        for (int i = 0; i != l; i++) {
            params.h[i].random_generator(get_random_bytes);
        }
    }

    void keygen(SecretKey& sk, const Params& params, const MasterKey& msk, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX rx;
        Scalar r;
        G1 temp;
        random_zpstar(rx, r, get_random_bytes);
        sk.a0.copy(params.g3);
        int j = 0; /* Index for writing to qualified.b */
        int k = 0; /* Index for reading from attrs.attrs */
        for (int i = 0; i != params.l; i++) {
            if (k != attrs.length && attrs.attrs[k].idx == i) {
                if (!attrs.attrs[k].omitFromKeys) {
                    temp.multiply(params.h[i], attrs.attrs[k].id);
                    sk.a0.add(sk.a0, temp);
                }
                k++;
            } else if (!attrs.omitAllFromKeysUnlessPresent) {
                sk.b[j].idx = i;
                sk.b[j].hexp.multiply(params.h[i], r);
                j++;
            }
        }
        sk.l = j;
        sk.signatures = params.signatures;
        if (sk.signatures) {
            sk.bsig.multiply(params.hsig, r);
        } else {
            sk.bsig.copy(G1::zero);
        }
        sk.a0.multiply(sk.a0, r);

        sk.a0.add(sk.a0, msk.g2alpha);
        sk.a1.multiply_frobenius(params.g, rx);
    }

    void qualifykey(SecretKey& qualified, const Params& params, const SecretKey& sk, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX tx;
        Scalar t;
        G1 temp;
        G1 product;
        random_zpstar(tx, t, get_random_bytes);
        product.copy(params.g3);
        qualified.a0.copy(sk.a0);
        int j = 0; /* Index for writing to qualified.b */
        int k = 0; /* Index for reading from attrs.attrs */
        int x = 0; /* Index for reading from sk.b */
        for (int i = 0; i != params.l; i++) {
            if (k != attrs.length && attrs.attrs[k].idx == i) {
                if (!attrs.attrs[k].omitFromKeys) {
                    temp.multiply(params.h[i], attrs.attrs[k].id);
                    product.add(product, temp);
                    if (x != sk.l && sk.b[x].idx == i) {
                        temp.multiply(sk.b[x].hexp, attrs.attrs[k].id);
                        qualified.a0.add(qualified.a0, temp);
                        x++;
                    }
                }
                k++;
            } else if (x != sk.l && sk.b[x].idx == i) {
                if (!attrs.omitAllFromKeysUnlessPresent) {
                    qualified.b[j].idx = i;
                    qualified.b[j].hexp.multiply(params.h[i], t);
                    qualified.b[j].hexp.add(qualified.b[j].hexp, sk.b[x].hexp);
                    j++;
                }
                x++;
            }
            /*
             * We could hit neither case if slot i is "hidden" (the b element
             * corresponding to it is not provided, so it can't be filled in).
             */
        }
        qualified.l = j;
        qualified.signatures = sk.signatures;
        if (qualified.signatures) {
            qualified.bsig.multiply(params.hsig, t);
            qualified.bsig.add(qualified.bsig, sk.bsig);
        } else {
            qualified.bsig.copy(G1::zero);
        }
        product.multiply(product, t);
        qualified.a0.add(qualified.a0, product);
        qualified.a1.multiply_frobenius(params.g, tx);
        qualified.a1.add(qualified.a1, sk.a1);
    }

    void nondelegable_keygen(SecretKey& sk, const Params& params, const MasterKey& msk, const AttributeList& attrs) {
        G1 temp;
        sk.a0.copy(params.g3);
        int j = 0; /* Index for writing to qualified.b */
        int k = 0; /* Index for reading from attrs.attrs */
        for (int i = 0; i != params.l; i++) {
            if (k != attrs.length && !attrs.attrs[k].omitFromKeys && attrs.attrs[k].idx == i) {
                temp.multiply(params.h[i], attrs.attrs[k].id);
                sk.a0.add(sk.a0, temp);
                k++;
            } else if (!attrs.omitAllFromKeysUnlessPresent) {
                sk.b[j].idx = i;
                sk.b[j].hexp.copy(params.h[i]);
                j++;
            }
        }
        sk.l = j;
        sk.signatures = params.signatures;
        if (sk.signatures) {
            sk.bsig.copy(params.hsig);
        } else {
            sk.bsig.copy(G1::zero);
        }
        sk.a0.add(sk.a0, msk.g2alpha);
        sk.a1.copy(params.g);
    }

    void nondelegable_qualifykey(SecretKey& qualified, const Params& params, const SecretKey& sk, const AttributeList& attrs) {
        G1 temp;
        qualified.a0.copy(sk.a0);
        int j = 0; /* Index for writing to qualified.b */
        int k = 0; /* Index for reading from attrs.attrs */
        int x = 0; /* Index for reading from sk.b */
        for (int i = 0; x != sk.l && i != params.l; i++) {
            if (k != attrs.length && attrs.attrs[k].idx == i) {
                if (sk.b[x].idx == i && !attrs.attrs[k].omitFromKeys) {
                    temp.multiply(sk.b[x].hexp, attrs.attrs[k].id);
                    qualified.a0.add(qualified.a0, temp);
                    x++;
                }
                k++;
            } else if (sk.b[x].idx == i) {
                if (!attrs.omitAllFromKeysUnlessPresent) {
                    qualified.b[j].idx = i;
                    qualified.b[j].hexp.copy(sk.b[x].hexp);
                    j++;
                }
                x++;
            }
            /*
             * We could hit neither case if slot i is "hidden" (the b element
             * corresponding to it is not provided, so it can't be filled in).
             */
        }
        qualified.l = j;
        qualified.signatures = sk.signatures;
        if (qualified.signatures) {
            qualified.bsig.copy(sk.bsig);
        } else {
            qualified.bsig.copy(G1::zero);
        }
        qualified.a1.copy(sk.a1);
    }

    void adjust_nondelegable(SecretKey& sk, const SecretKey& parent, const AttributeList& from, const AttributeList& to) {
        G1 temp;
        Scalar diff;

        int j = 0;
        int k = 0;
        int x = 0;
        for (int i = 0; i != parent.l; i++) {
            int idx = parent.b[i].idx;
            while (j != from.length && from.attrs[j].idx < idx && !from.attrs[j].omitFromKeys) {
                j++;
            }
            while (k != to.length && to.attrs[k].idx < idx && !to.attrs[k].omitFromKeys) {
                k++;
            }

            bool sub_from = (j != from.length && from.attrs[j].idx == idx);
            bool add_to = (k != to.length && to.attrs[k].idx == idx);

            if (j != from.length || k != to.length) {
                if (sub_from && add_to) {
                    if (!ID::equal(from.attrs[j].id, to.attrs[k].id)) {
                        if (diff.subtract(to.attrs[k].id, from.attrs[j].id)) {
                            diff.add(diff, group_order);
                        }
                        temp.multiply(parent.b[i].hexp, diff);
                        sk.a0.add(sk.a0, temp);
                    }
                } else if (sub_from) {
                    diff.subtract(group_order, from.attrs[j].id);
                    temp.multiply(parent.b[i].hexp, diff);
                    sk.a0.add(sk.a0, temp);
                } else if (add_to) {
                    temp.multiply(parent.b[i].hexp, to.attrs[k].id);
                    sk.a0.add(sk.a0, temp);
                }
            }

            if (!add_to) {
                sk.b[x].idx = parent.b[i].idx;
                sk.b[x].hexp.copy(parent.b[i].hexp);
                x++;
            }
        }

        sk.l = x;
    }

    void precompute(Precomputed& precomputed, const Params& params, const AttributeList& attrs) {
        G1 temp;
        precomputed.prodexp.copy(params.g3);
        for (int i = 0; i != attrs.length; i++) {
            const Attribute& attr = attrs.attrs[i];
            temp.multiply(params.h[attr.idx], attr.id);
            precomputed.prodexp.add(precomputed.prodexp, temp);
        }
    }

    void adjust_precomputed(Precomputed& precomputed, const Params& params, const AttributeList& from, const AttributeList& to) {
        G1 temp;
        Scalar diff;

        int i = 0;
        int j = 0;
        while (i != from.length && j != to.length) {
            const Attribute& from_attr = from.attrs[i];
            const Attribute& to_attr = to.attrs[j];
            if (from_attr.idx == to_attr.idx) {
                if (!ID::equal(from_attr.id, to_attr.id)) {
                    if (diff.subtract(to_attr.id, from_attr.id)) {
                        diff.add(diff, group_order); // b-a
                    }
                    temp.multiply(params.h[to_attr.idx], diff);
                    precomputed.prodexp.add(precomputed.prodexp, temp);
                }
                i++;
                j++;
            } else if (from_attr.idx < to_attr.idx) {
                diff.subtract(group_order, from_attr.id);
                temp.multiply(params.h[from_attr.idx], diff);
                precomputed.prodexp.add(precomputed.prodexp, temp); // -ai
                i++;
            } else {
                temp.multiply(params.h[to_attr.idx], to_attr.id); // bi
                precomputed.prodexp.add(precomputed.prodexp, temp);
                j++;
            }
        }
        while (i != from.length) {
            const Attribute& from_attr = from.attrs[i];
            diff.subtract(group_order, from_attr.id);
            temp.multiply(params.h[from_attr.idx], diff);
            precomputed.prodexp.add(precomputed.prodexp, temp);
            i++;
        }
        while (j != to.length) {
            const Attribute& to_attr = to.attrs[j];
            temp.multiply(params.h[to_attr.idx], to_attr.id);
            precomputed.prodexp.add(precomputed.prodexp, temp);
            j++;
        }
    }

    void resamplekey(SecretKey& resampled, const Params& params, const Precomputed& precomputed, const SecretKey& sk, bool supportFurtherQualification, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX tx;
        Scalar t;
        G1 temp;
        G2 temp2;
        random_zpstar(tx, t, get_random_bytes);

        temp.multiply(precomputed.prodexp, t);
        resampled.a0.add(sk.a0, temp);

        temp2.multiply_frobenius(params.g, tx);
        resampled.a1.add(sk.a1, temp2);

        resampled.signatures = sk.signatures;
        if (resampled.signatures) {
            temp.multiply(params.hsig, t);
            resampled.bsig.add(sk.bsig, temp);
        } else {
            resampled.bsig.copy(G1::zero);
        }

        if (supportFurtherQualification) {
            for (int i = 0; i != sk.l; i++) {
                temp.multiply(params.h[sk.b[i].idx], t);
                resampled.b[i].hexp.add(sk.b[i].hexp, temp);
                resampled.b[i].idx = sk.b[i].idx;
            }
            resampled.l = sk.l;
        } else {
            resampled.l = 0;
        }
    }

    void encrypt(Ciphertext& ciphertext, const GT& message, const Params& params, const AttributeList& attrs, void (*get_random_bytes)(void*, size_t)) {
        Precomputed precomputed;
        precompute(precomputed, params, attrs);
        encrypt_precomputed(ciphertext, message, params, precomputed, get_random_bytes);
    }

    void encrypt_precomputed(Ciphertext& ciphertext, const GT& message, const Params& params, const Precomputed& precomputed, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX sx;
        Scalar s;
        random_zpstar(sx, s, get_random_bytes);

        ciphertext.a.exponentiate_gt(params.pairing, sx);
        ciphertext.a.multiply(ciphertext.a, message);
        ciphertext.b.multiply_frobenius(params.g, sx);
        ciphertext.c.multiply(precomputed.prodexp, s);
    }

    void decrypt(GT& message, const Ciphertext& ciphertext, const SecretKey& sk) {
        GT denominator;
        G1Affine caffine;
        G2Affine a1affine;
        G1Affine a0affine;
        G2Affine baffine;
        caffine.from_projective(ciphertext.c);
        a1affine.from_projective(sk.a1);
        a0affine.from_projective(sk.a0);
        baffine.from_projective(ciphertext.b);

        a0affine.negate(a0affine);
        bls12_381::AffinePair pairs[2];
        pairs[0].g1 = &caffine;
        pairs[0].g2 = &a1affine;
        pairs[1].g1 = &a0affine;
        pairs[1].g2 = &baffine;
        bls12_381::pairing_product(message, pairs, 2, nullptr, 0);
        message.multiply(message, ciphertext.a);
    }

    void decrypt_master(GT& message, const Ciphertext& ciphertext, const MasterKey& msk) {
        G1Affine g2alphaaffine;
        G2Affine baffine;
        g2alphaaffine.from_projective(msk.g2alpha);
        baffine.from_projective(ciphertext.b);

        g2alphaaffine.negate(g2alphaaffine);
        bls12_381::pairing(message, g2alphaaffine, baffine);
        message.multiply(message, ciphertext.a);
    }

    void sign(Signature& signature, const Params& params, const SecretKey& sk, const AttributeList* attrs, const Scalar& message, void (*get_random_bytes)(void*, size_t)) {
        Precomputed precomputed;
        precompute(precomputed, params, *attrs);
        sign_precomputed(signature, params, sk, attrs, precomputed, message, get_random_bytes);
    }

    void sign_precomputed(Signature& signature, const Params& params, const SecretKey& sk, const AttributeList* attrs, const Precomputed& precomputed, const Scalar& message, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX sx;
        Scalar s;
        G1 prodexp;
        random_zpstar(sx, s, get_random_bytes);

        signature.a0.multiply(sk.bsig, message);
        prodexp.multiply(params.hsig, message);
        signature.a0.add(signature.a0, sk.a0);
        prodexp.add(prodexp, precomputed.prodexp);
        signature.a1.multiply_frobenius(params.g, sx);
        prodexp.multiply(prodexp, s);
        signature.a0.add(signature.a0, prodexp);
        signature.a1.add(signature.a1, sk.a1);

        if (attrs != nullptr) {
            int k = 0;
            for (int i = 0; i != sk.l; i++) {
                while (k != attrs->length && attrs->attrs[k].idx < sk.b[i].idx) {
                    k++;
                }
                if (k == attrs->length) {
                    return;
                }
                if (sk.b[i].idx == attrs->attrs[k].idx) {
                    prodexp.multiply(sk.b[i].hexp, attrs->attrs[k].id);
                    signature.a0.add(signature.a0, prodexp);
                    k++;
                }
            }
        }
    }

    bool verify(const Params& params, const AttributeList& attrs, const Signature& signature, const Scalar& message) {
        Precomputed precomputed;
        precompute(precomputed, params, attrs);
        return verify_precomputed(params, precomputed, signature, message);
    }

    bool verify_precomputed(const Params& params, const Precomputed& precomputed, const Signature& signature, const Scalar& message) {
        G1Affine a0affine;
        G2Affine gaffine;
        G1Affine prodexpaffine;
        G2Affine a1affine;

        {
            G1 prodexp;
            prodexp.multiply(params.hsig, message);
            prodexp.add(prodexp, precomputed.prodexp);
            a0affine.from_projective(signature.a0);
            gaffine.from_projective(params.g);
            prodexpaffine.from_projective(prodexp);
            a1affine.from_projective(signature.a1);
        }

        /* Compute e(a0affine, gaffine) / e(prodexpaffine, a1affine). */
        GT ratio;
        prodexpaffine.negate(prodexpaffine);
        bls12_381::AffinePair pairs[2];
        pairs[0].g1 = &a0affine;
        pairs[0].g2 = &gaffine;
        pairs[1].g1 = &prodexpaffine;
        pairs[1].g2 = &a1affine;
        bls12_381::pairing_product(ratio, pairs, 2, nullptr, 0);

        return GT::equal(ratio, params.pairing);
    }

    std::vector<std::vector<uint8_t>> KeyDerivationTree::generate_keys(const uint8_t* seed, size_t seed_len, size_t num_keys, size_t key_len) {
        std::vector<std::vector<uint8_t>> keys;
        std::vector<uint8_t> current(seed, seed + seed_len);

        for (size_t i = 0; i < num_keys; ++i) {
            std::vector<uint8_t> key(key_len);
            for (size_t j = 0; j < key_len; ++j) {
                key[j] = current[j % current.size()] ^ static_cast<uint8_t>(i + j);
            }
            keys.push_back(key);
            current = key;
        }
        return keys;
    }

    size_t SymmetricCipher::encrypt(uint8_t* out, const uint8_t* data, size_t len, const uint8_t* key, size_t key_len) {
        if (key_len != 32) {
            throw std::runtime_error("SymmetricCipher::encrypt: key must be 32 bytes for AES-256");
        }

        EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
        if (!ctx) throw std::runtime_error("Failed to create EVP_CIPHER_CTX");

        uint8_t iv[16];
        if (!RAND_bytes(iv, sizeof(iv))) throw std::runtime_error("Failed to generate IV");
        memcpy(out, iv, sizeof(iv));

        int out_len1 = 0, out_len2 = 0;
        if (!EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), nullptr, key, iv))
            throw std::runtime_error("EncryptInit failed");

        if (!EVP_EncryptUpdate(ctx, out + 16, &out_len1, data, len))
            throw std::runtime_error("EncryptUpdate failed");

        if (!EVP_EncryptFinal_ex(ctx, out + 16 + out_len1, &out_len2))
            throw std::runtime_error("EncryptFinal failed");

        EVP_CIPHER_CTX_free(ctx);
        return 16 + out_len1 + out_len2;
    }

    size_t SymmetricCipher::decrypt(uint8_t* out, const uint8_t* ciphertext, size_t len, const uint8_t* key, size_t key_len) {
        if (key_len != 32 || len < 16)
            throw std::runtime_error("SymmetricCipher::decrypt: invalid input");

        const uint8_t* iv = ciphertext;
        const uint8_t* enc = ciphertext + 16;
        int enc_len = len - 16;

        EVP_CIPHER_CTX* ctx = EVP_CIPHER_CTX_new();
        if (!ctx) throw std::runtime_error("Failed to create EVP_CIPHER_CTX");

        int out_len1 = 0, out_len2 = 0;
        if (!EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), nullptr, key, iv))
            throw std::runtime_error("DecryptInit failed");

        if (!EVP_DecryptUpdate(ctx, out, &out_len1, enc, enc_len))
            throw std::runtime_error("DecryptUpdate failed");

        if (!EVP_DecryptFinal_ex(ctx, out + out_len1, &out_len2))
            throw std::runtime_error("DecryptFinal failed");

        EVP_CIPHER_CTX_free(ctx);
        return out_len1 + out_len2;
    }

    static void write_uint32(std::vector<uint8_t>& out, uint32_t val) {
        for (int i = 3; i >= 0; --i)
            out.push_back((val >> (i * 8)) & 0xFF);
    }

    static uint32_t read_uint32(const uint8_t*& ptr) {
        uint32_t val = 0;
        for (int i = 0; i < 4; ++i)
            val = (val << 8) | *ptr++;
        return val;
    }

    static void write_uint64(std::vector<uint8_t>& out, uint64_t val) {
        for (int i = 7; i >= 0; --i)
            out.push_back((val >> (i * 8)) & 0xFF);
    }

    static uint64_t read_uint64(const uint8_t*& ptr) {
        uint64_t val = 0;
        for (int i = 0; i < 8; ++i)
            val = (val << 8) | *ptr++;
        return val;
    }

    static void write_bytes(std::vector<uint8_t>& out, const std::vector<uint8_t>& data) {
        write_uint32(out, data.size());
        out.insert(out.end(), data.begin(), data.end());
    }

    static std::vector<uint8_t> read_bytes(const uint8_t*& ptr) {
        uint32_t len = read_uint32(ptr);
        return std::vector<uint8_t>(ptr, ptr + len);
        ptr += len;
    }

    std::vector<uint8_t> TokenEncryptor::serialize(const Token& t) {
        std::vector<uint8_t> out;

        auto append_u32 = [&](uint32_t val) {
            for (int i = 3; i >= 0; --i)
                out.push_back((val >> (i * 8)) & 0xFF);
        };
        auto append_u64 = [&](uint64_t val) {
            for (int i = 7; i >= 0; --i)
                out.push_back((val >> (i * 8)) & 0xFF);
        };
        auto append_vec = [&](const std::vector<uint8_t>& v) {
            append_u32(v.size());
            out.insert(out.end(), v.begin(), v.end());
        };
        auto append_str = [&](const std::string& s) {
            append_u32(s.size());
            out.insert(out.end(), s.begin(), s.end());
        };

        append_u32(t.id);
        append_str(t.prg);
        append_str(t.uid);
        append_u32(t.height);
        append_u32(t.offset);
        append_vec(t.seed);
        append_u64(t.gen_time);
        append_u64(t.exp_time);
        append_vec(t.auth);
        return out;
    }

    Token TokenEncryptor::deserialize(const std::vector<uint8_t>& buffer) {
        Token t;
        const uint8_t* ptr = buffer.data();
        const uint8_t* end = buffer.data() + buffer.size();

        auto read_u32 = [&]() -> uint32_t {
            if (end - ptr < 4) throw std::runtime_error("deserialize: not enough bytes for u32");
            uint32_t val = 0;
            for (int i = 0; i < 4; ++i) val = (val << 8) | *ptr++;
            return val;
        };
        auto read_u64 = [&]() -> uint64_t {
            if (end - ptr < 8) throw std::runtime_error("deserialize: not enough bytes for u64");
            uint64_t val = 0;
            for (int i = 0; i < 8; ++i) val = (val << 8) | *ptr++;
            return val;
        };
        auto read_bytes = [&](std::vector<uint8_t>& v) {
            uint32_t len = read_u32();
            if (end - ptr < len) throw std::runtime_error("deserialize: byte vector out of bounds");
            v.assign(ptr, ptr + len);
            ptr += len;
        };
        auto read_string = [&](std::string& s) {
            uint32_t len = read_u32();
            if (end - ptr < len) throw std::runtime_error("deserialize: string out of bounds");
            s.assign((const char*)ptr, len);
            ptr += len;
        };

        t.id = read_u32();
        read_string(t.prg);
        read_string(t.uid);
        t.height = read_u32();
        t.offset = read_u32();
        read_bytes(t.seed);
        t.gen_time = read_u64();
        t.exp_time = read_u64();
        read_bytes(t.auth);
        uint32_t prg_len = read_u32();
        printf("[debug] deserialize: prg_len = %u, remaining = %ld\n", prg_len, end - ptr);
        return t;
    }


    static void hash_to_gt(const uint8_t* key, size_t key_len, GT& out_gt) {
        uint8_t buffer[384] = {0};
        uint8_t hash[SHA256_DIGEST_LENGTH];

        SHA256(key, key_len, hash);
        memcpy(buffer, hash, std::min(sizeof(buffer), sizeof(hash)));

        out_gt.read_big_endian(buffer);
    }

    void default_hash_gt_to_key(uint8_t* out_key, size_t key_len, const GT& gt) {
        uint8_t buffer[384];
        gt.write_big_endian(buffer);

        uint8_t hash[SHA256_DIGEST_LENGTH];
        SHA256(buffer, sizeof(buffer), hash);

        memcpy(out_key, hash, key_len);
    }

    void TokenEncryptor::encrypt_token(Ciphertext& out_ct,
                                       std::vector<uint8_t>& out_sym_cipher,
                                       const Token& token,
                                       const Params& p,
                                       const AttributeList& attrs,
                                       void (*random_bytes)(void*, size_t)) {
        uint8_t sym_key[32];
        random_bytes(sym_key, sizeof(sym_key));

        std::vector<uint8_t> plaintext = serialize(token);
        out_sym_cipher.resize(plaintext.size());
        SymmetricCipher::encrypt(out_sym_cipher.data(), plaintext.data(), plaintext.size(), sym_key, sizeof(sym_key));

        GT gt_key;
        hash_to_gt(sym_key, sizeof(sym_key), gt_key);  // 对称 hash 生成 GT

        encrypt(out_ct, gt_key, p, attrs, random_bytes);
    }

    Token TokenEncryptor::decrypt_token(const Ciphertext& ct,
                                        const std::vector<uint8_t>& sym_cipher,
                                        const SecretKey& sk) {

        GT gt_key;
        decrypt(gt_key, ct, sk);

        uint8_t sym_key[32];
        default_hash_gt_to_key(sym_key, sizeof(sym_key), gt_key);

        std::vector<uint8_t> decrypted(sym_cipher.size());
        SymmetricCipher::decrypt(decrypted.data(), sym_cipher.data(), sym_cipher.size(), sym_key, sizeof(sym_key));

        return deserialize(decrypted);
    }


}
