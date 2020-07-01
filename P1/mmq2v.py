import numpy as np

u = [1, 1, 1, 1, 1]
v = [2, 4, 6, 8, 10]
b1 = np.array([np.log(51.34), np.log(52.72), np.log(40.6), np.log(27.769), np.log(17.84)])
b2 = np.array([np.log(2), np.log(4), np.log(6), np.log(8), np.log(10)])
b=b1-b2

u = np.array(u)
v = np.array(v)
b = np.array(b)

uu = u.dot(u)
uv = u.dot(v)
vv = v.dot(v)
ub = u.dot(b)
vb = v.dot(b)

dire = np.array([[uu, uv], [uv, vv]])
esq = np.array([ub, vb])
print('Resolvendo:')
print(dire, '=', esq)
res = np.linalg.solve(dire, esq)
print('Res', res)
