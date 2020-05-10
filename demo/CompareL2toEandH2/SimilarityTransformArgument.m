%% System
G = rss(3,2,4);
G.D = 0;
[A,B,C,D] = ssdata(G);

%%
Gfh = evalt(tvss(G),0:0.1:1);
tvnorm(Gfh,2)

%%
[U,S,V] = svd(G.C,0);
G1 = ss(V*A*V',V*B,C*V',0);
G1fh = evalt(tvss(G1),0:0.1:1);
tvnorm(G1fh,2)