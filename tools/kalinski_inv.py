
import math

def random_prime(bound):
    warnings.simplefilter('ignore')
    chck = False
    while chck == False:
        p = random.randrange(bound>>1, bound)
        chck = pyprimes.isprime(p)
    warnings.simplefilter('default')
    return p

def fieldElement(a):
	f0 = a % pow(2,52)
	a = a / pow(2,52)
	f1 = a % pow(2,52)
	a = a / pow(2,52)
	f2 = a % pow(2,52)
	a = a / pow(2,52)
	f3 = a % pow(2,52)
	a = a / pow(2,52)
	f4 = a % pow(2,52)
	arr = [f0, f1, f2, f3, f4]
	return '[' + ', '.join([str(x) for x in arr]) + ']'


def AlmostMontInv(x, p):
    u = p
    v = x
    r = 0
    s = 1
    k = 0
    while (v>0):
        if (u%2 == 0):
            u = u >> 1
            s = s << 1
        elif (v%2 == 0):
            v = v >> 1
            r = r << 1
        elif u > v:
            u = (u-v) >> 1
            r = r+s
            s = s << 1
        else:
            v = (v-u) >> 1
            s = r+s
            r = r << 1
        k += 1
        #print('Values on iteration: ' + str(k) + ':')
        #print('r = ' + str(fieldElement(r)))
        #print('s = ' + str(fieldElement(s)))
        #print('v = ' + str(fieldElement(v)))
        #print('u = ' + str(fieldElement(u)))
        
    if r>p:
        r = r-p

    return p-r, k


def phase2(r, k):
    for i in range(0, k-253):
        if (r % 2 == 0):
            r = (r >> 1) % p
        else:
            r = (r+p) >> 1
        #print('r on iter: ' + str(i) + '= ' + str(r))
        print('iter: ' + str(i) + ': FieldElement(' + fieldElement(r) + ')')
    return r

# random.seed(11) # for debug

p = 7237005577332262213973186563042994240857116359379907606001950938285454250989

# Input to change value to get the AlmostMontInv.
a = 904625697166532776746648320197686575422163851717637391703244652875051672039

ainv, k = AlmostMontInv(a, p)
print('r:FieldElement' + str(fieldElement(ainv)))
print('r:FieldElement' + str(ainv))
print('k:' + str(k))
print('p:' + str(fieldElement(p)))
ainv = phase2(ainv, k)

print(ainv, k)
print('ainv:FieldElement' + str(fieldElement(ainv)))
