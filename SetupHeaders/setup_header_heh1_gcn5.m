% Setup header/info for double KO heh1D gcn5D nuclear fluctuations.

namegroup1 = {'wt','gcn5','heh1_gcn5','heh1'};
nameIntercleave = reshape([namegroup1],1,1*length(namegroup1));
nameAll = [namegroup1];
nameAll2 = {'WT','\itgcn5\Delta','\itheh1\Deltagcn5\Delta',...
    '\itheh1\Delta'};
colorAll = [0 0 0; .91,.41,.17; .52,.34,.14; 1 1 0];
numM=length(namegroup1);
[ points2,faces2,psi,theta] = setup_img( );
p2um=0.16;
f2s=2.5;
[points,faces,edges,neighbors] = TriSphere(3);