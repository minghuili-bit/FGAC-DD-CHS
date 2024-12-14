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

using namespace std;

namespace embedded_pairing::wkdibe {

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

    void setup(Params& params, MasterKey& msk, int l, bool signatures, void (*get_random_bytes)(void*, size_t)) {
        bls12_381::PowersOfX alphax;
        Scalar alpha;
        random_zpstar(alphax, alpha, get_random_bytes); // 使用alphax函数生成一个随机的标量alpha
        params.g.random_generator(get_random_bytes);  // 生成一个随机基数g
        params.g1.multiply_frobenius(params.g, alphax);
        params.g2.random_generator(get_random_bytes);

        std::cout << "Value of alpha: " << static_cast<int>(alpha) << std::endl;

        // Vector to store the points
        vector<pair<int, int> > points(4);

        // Sharing of secret Code in N parts
        secret_sharing(static_cast<int>(alpha), points, 4, 2);

        cout << "Secret is divided to " << 4
             << " Parts - " << endl;

        for (int i = 0; i < 4; ++i) {
            cout << points[i].first << " "
                 << points[i].second << endl;
            msk.points.push_back(points[i]);
        }

        for (int i = 0; i < 4; i++) {
            auto secret = Scalar{.std_words = {(uint32_t)points[i].second}};
            std::cout << "Value of alpha: " << static_cast<int>(secret) << std::endl;
            msk.g2alpha1[i].multiply(params.g2, secret);
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
}
