
def compositions(n,k):
  '''
  Creates a generator for all compositions of the integer n as k summands.
  E.g. compositions(4,3) gives (0,0,4), (0,3,1), ..., (4,0,0).
  Different orderings are counted as different compositions. Code taken
  from https://dandrake.livejournal.com/83095.html

  :param n: Positive integer to be composed
  :param k: Number of summands in the composition
  :return: Generator giving all compositions as tuples
  '''
  if n < 0 or k < 0:
    return
  elif k == 0:
    # the empty sum, by convention, is zero, so only return something if
    # n is zero
    if n == 0:
      yield []
    return
  elif k == 1:
    yield [n]
    return
  else:
    for i in range(0,n+1):
      for comp in compositions(n-i,k-1):
        yield [i] + comp

def test_compositions():
  print([k for k in compositions(4, 3)])
  print(len([k for k in compositions(4, 3)]))


if __name__ == '__main__':
    test_compositions()
