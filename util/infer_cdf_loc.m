function p=infer_cdf_loc(Y,x1)

[f,x]=ecdf(Y);

x_dis=abs(x-x1);

x_inCdf=x(find(x_dis==min(x_dis)));
f_inCdf=f(find(x_dis==min(x_dis)));
p=mean([x_inCdf,f_inCdf],1);