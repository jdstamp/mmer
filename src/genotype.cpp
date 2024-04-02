/*
 * Substantial parts of this file were published under MIT license by other
 * authors. Find the original license notice below.
 */

/* MIT License
 *
 * Copyright (c) 2024 sriramlab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "genotype.h"

using namespace std;

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value) {
  return stream.read(reinterpret_cast<char *>(&value), sizeof(T));
}

double genotype::get_col_mean(int snpindex) const {
  double temp = columnmeans[snpindex];
  return temp;
}

double genotype::get_col_std(int snpindex) const {
  double p_i = get_col_mean(snpindex);
  if (p_i == 0 || p_i == 2)
    return 1.0;
  double temp = sqrt(p_i * (1 - (0.5 * p_i)));
  return temp;
}
