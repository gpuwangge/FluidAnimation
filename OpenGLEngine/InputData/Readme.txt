Code format:
horizontal direction: width, x, i
vertical direction: height, y, j

v[y][x]
v[j][i]


Format of Config.txt
=========================
#1. Resolution of the grid
0: 4x8
1: 6x12
2: 32x64
3: 64x128
4: 128x256
5: 256x512
6: 25x50

#2. Max step
Set to the value according to the input velocity field
For airfoil-velocity.dat, there is 91 steps
If it is zero, use generated velocity field instead


#3. Should it pause after each frame
0: No
1: Yes

#4. Down sample
Can not be zero

#5. Show unit vector
0: No
1: Yes


=========================
Format of Parameters.txt
#1. Integration T
#2. Integration dt
#3. Epsilon: S = 1/(t^2+epsilon)
#4. Threshold: if(lambda_min > threshold) S = 1
#5. Epsilon: t = s^2 *ee^T + epsilon*I
#6. Weight: B = L'L + weight*D
#7. Smooth h: (M - h * L) * S^(t+1) = M * S^t
#8. Smooth iteration n
#9. Epsilon: if (FTLE < epsilon) S = 1; else S = 5
(for generated test, this is 0.3; for airfoil, this is 0.6667?)



=========================

for generated test: use T=15, dt = 0.05, maxstep = 0
for airfoil: use T = -4, dt = 0.1, maxstep = 91