import numpy as np

u = [1, 1, 1, 1, 1]
v = [2, 4, 6, 8, 10]
w = []
b1 = np.array([np.log(51.34), np.log(52.72), np.log(40.6), np.log(27.769), np.log(17.84)])
b2 = np.array([np.log(2), np.log(4), np.log(6), np.log(8), np.log(10)])
b=b1-b2

u = np.array(u)
v = np.array(v)
w = np.array(w)
b = np.array(b)

uu = u.dot(u)
uv = u.dot(v)
uw = u.dot(w)
vv = v.dot(v)
vw = v.dot(w)
ww = v.dot(w)
ub = u.dot(b)
vb = v.dot(b)
wb = w.dot(b)

dire = np.array([[uu, uv, uw], [uv, vv, vw], [uw, vw, ww]])
esq = np.array([ub, vb, wb])
print('Resolvendo:')
print(dire, '=', esq)
res = np.linalg.solve(dire, esq)
print('Res', res)
