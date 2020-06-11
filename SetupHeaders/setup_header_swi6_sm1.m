namegroup1={'wt','swi6','swi6-sm1'};
% namegroup2={'wt_MBC','swi6_MBC','epe1'};
nameIntercleave=reshape([namegroup1],1,...
    1*length(namegroup1));
% nameIntercleave=reshape([namegroup1;namegroup2],1,...
%     2*length(namegroup1));
nameAll=namegroup1;
nameAll2=cellfun(@(x)strrep(x,'_',' '),nameAll,'UniformOutput',0);
colorAll=[0,0,0;1,0,0;.39, .2, .6;0, 0, 0;1, 0, 0;0 0 1];
numM=length(namegroup1);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);