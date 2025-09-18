@namespace "random"

#' Generate a random integer between 1 and `n`, inclusive.
#'
#' @param n Maximum.
#' @return Random integer between 1 and `n`, inclusive. If `n` is not positive,
#'   the function will return 0.
function unif(n) {
	if (n <= 0) {
		return 0
	}

	return 1 + int(int(n) * rand())
}

#' Sample `m` integers between 1 and `n`, inclusive.
#'
#' @param m Number of samples to take. If `m >= n`, then `n` samples will be
#'   taken. If `m` or `n` is non-positive, no samples will be taken.
#' @param n Maximum.
#' @param arr Samples will be assigned as values in `arr`. The order of the
#'   assignments is arbitrary and `arr` will be cleared before assignment.
#' @param [replace] If 0 or "" (the default) samples are taken without
#'   replacement. Any other value will cause samples to be taken with
#'   replacement.
#' @return Number of samples taken.
function sample(m, n, arr, replace) {
	delete arr

	if (m <= 0 || n <= 0) {
		return 0
	}

	m = int(m)
	n = int(n)
	if (m >= n) {
		m = n
	}

	if (replace) {
		return _sample_replace(m, n, arr)
	}

	return _sample_floyd(m, n, arr)
}

#' Sample `m` indicies from an array.
#'
#' Often, it is useful to sample from an array in which the indicies are not
#' integers 1 to n. In order to sample from such an array with randomly
#' generated integers, all the indicies in the array must be copied to a
#' temporary array indexed with integers and then sampled. A more memory
#' efficient approach is to use reservoir sampling. First `m` indicies are
#' copied from `src` to `dest` (assigned to values of `dest`). Then for each
#' ith index in `src` beyond the `m`th, it replaces an entry in `dest` with '
#' probability m / i.
#'
#' @param m Number of samples to take.
#' @param src Array from which to take samples. The indicies are sampled.
#' @param dest Array that will contain the sampled indices. `dest` will be
#'   indexed from 1 to `m` or the length of `src`, whichever is less.
#' @return Number of samples taken.
function sample_array(m, src, dest,    x, i, j) {
	delete dest
	m = int(m)
	for (x in src) {
		if (++i <= m) {
			dest[i] = x
			continue
		}

		j = unif(i)
		if (j <= m) {
			dest[j] = x
		}
	}

	return i < m ? i : m
}

function _sample_replace(m, n, arr,    i) {
	for (i = 1; i <= m; ++i) {
		arr[i] = unif(n)
	}

	return m
}

# Robert Floyd's algorithm for random sampling.
# This implementation is copied almost verbatim from Jon Bentley's
# "programming pearls":
# "Communications of the ACM" September 1987, Volume 30, Number 9
# (PROGRAM F2)
# Read here "https://fermatslibrary.com/s/a-sample-of-brilliance"
#
# The algorithm starts with an empty set of selected values and a set of
# candidate values with a cardinality equal to n - m. At each iteration, a new
# value is added to the candidates set and then a random value from the set is
# added to the selected set. This continues until the selected set has m
# values. There is likely some complicated math needed prove that sampling is
# truly random, so I am going to trust Mr. Floyd on this.
#
# The algorithm is similar to the Fisher-Yates shuffle in that it can be viewed
# in terms of generating a permutation.
#
# +---------------------------------------------------------------------------+
# |    "rejected" values     |   "selected" values  | |     unseen values     |
# +---------------------------------------------------------------------------+
#  ^                          ^                      ^                       ^
#  1                      n - m + 1                  j                       n
#
# At each iteration, the value at j is "swapped" with a value in rejected set
# or itself. Once a value is added to the selected set, it is never removed.
# The last m values become the sampled set.
#
# The "programming pearls" article mentions the original algorithm as being a
# recursive one. Selecting a random subset of 5 numbers between 1 and 10 is
# equivalent to selecting a random subset of 4 numbers between 1 and 9,
# generating a random number between 1 and 10, then adding 10 to the subset if
# the random number is 10 or one of the four numbers in the subset. 10 is added
# with probability 5 / 10 which is equal to (4 choose 9) / (5 choose 10).
#
# In the general case, the number of subsets of size m that can be made from a
# set of n numbers is (n choose m). The number of subsets that include a given
# number i is (n - 1 choose m - 1). Thus the probability of a number i being in
# a random subset of size m is (n - 1 choose m - 1) / (n choose m).
function _sample_floyd(m, n, arr,    i, j, k, s) {
	delete s
	for (j = n - m + 1; j <= n; ++j) {
		t = unif(j)
		if (t in s) {
			s[j] = 1
		} else {
			s[t] = 1
		}
	}

	k = 0
	for (i in s) {
		arr[++k] = i
	}

	return m
}
