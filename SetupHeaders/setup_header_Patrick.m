namegroup1={'wt'};
namegroup2={'heh1qd'};
namegroup3={'heh1D'};
% namegroup4={'lag1'};
% namegroup5={'pombe'};
nameIntercleave=reshape([namegroup1;namegroup2;namegroup3],1,3*length(namegroup1));
nameAll=[namegroup1,namegroup2,namegroup3];
nameAll2=cellfun(@(x)strrep(x,'_',' '),nameAll,'UniformOutput',0);
colorAll=[0, 0, 0;0, 0, 1;1, 0, 0;.13, .54, .13;.91, .41, .17;...
        .39, .2, .6];
numM=length(namegroup1);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);