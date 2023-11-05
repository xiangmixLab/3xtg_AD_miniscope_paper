% Jaccard index for cluster similarity

function [J,R,D,AR,RR,od,s]=new_cluster_similarity_index_031523(p1,p2)

% 1. consensus matrix
[cm,dm,cm_p]=consensus_mat_calculation([p1,p2]);

% 2. necessary statistics
a=sum(sum(cm==2)); % pairs always stay together in p1, p2
d=sum(sum(dm==2)); % pairs always not stay together in p1, p2
b=sum(sum((cm_p{1}==1).*(cm_p{2}==0))); % pairs stay togehter in p1 but not p2
c=sum(sum((cm_p{1}==0).*(cm_p{2}==1))); % pairs stay togehter in p2 but not p1

% Jaccard
J=a/(a+b+c);

% Rand
R=(a+d)/(a+b+c+d);

% DICE
D=2*a/(2*a+b+c);

% Russel-Rao
RR=a/(a+b+c+d);

% Adjusted Rand
%AR=2*(a*d-b*c)/((a+b)*(d+b)+(a+c)*(d+c));
%AR=((a+b+c+d)*(a+d)-((a+b)*(a+c)+(c+d)*(b+d)))/((a+b+c+d)^2-((a+b)*(a+c)+(c+d)*(b+d)));
E=((a+b)*(a+c)+(c+d)*(b+d))/(a+b+c+d)^2;
AR=(R-E)/(1-E);

% old
od=new_cluster_overlap_latest(p1,p2);
% statistic
s=[a,b,c,d];