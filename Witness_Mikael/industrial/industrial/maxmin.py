landmarks = []
landmarks.append(idx[0])
pt = tak[[landmarks[-1]],:]
N = 150
dists = zeros((N,tak.shape[0]))
distmin = spd.cdist(pt, tak)
dists[0,:] = distmin

for j in range(1,N):
    landmarks.append(distmin.argmax())
    pt = tak[[landmarks[-1]],:]
    dists[j,:] = spd.cdist(pt,tak)
    distmin = dists[:j+1,:].min(axis=0)
