from initializetracks import initializeTracks

num_azim = 16 #number of azimuthal angles desired
t = 0.05#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 2 #number of polar divisions; can be 2 or 3

tracks = initializeTracks(num_azim, t, w, h, n_p)
tracks.getTracks()