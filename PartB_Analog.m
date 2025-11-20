clc
%%Analog Communication Part B (FM)%%
n_max=[];
beta=0:0.2:40;
for j=beta

for i=1:40
  if(abs(besselj(40-i,j))>=0.01)
  break;
  end;
  endfor;
n_max(end+1)=40-i;

endfor
plot(beta,2.*n_max./beta)   % (normalized BW)= BW/delta_f = 2*n_max/beta
