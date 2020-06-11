% rootpath='.';
namegroup1={'wt'};
namegroup2={'wt MBC'};
% nameIntercleave=reshape(namegroup1,1,2*length(namegroup1));
nameIntercleave=namegroup1;
nameAll=namegroup1;
nameAll2=cellfun(@(x)strrep(x,'_',' '),nameAll,'UniformOutput',0);
colorAll=[255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162;255 0 0; 247 130 41; 255 255 0; 133 196 65; 17 106 54;...
    112 205 221; 58 83 164; 123 83 162]/255;
numM=length(namegroup1);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);