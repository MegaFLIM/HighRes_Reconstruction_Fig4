%%%%function for sliding window spatial binning.Window is 3 by 3 and is centerd at the
%%%%pixel. The input A is the matrix on which the binning is to be done. k is the size of the window. 
function Binned = SlideWinSum(A,k)
k = k-2;
row = size(A,1);
col = size(A,2);
A_horz = [zeros(row,1) A zeros(row,1)];
A_vert = [zeros(1,col+2); A_horz; zeros(1,col+2)];
Binned = zeros(row,col);
for i =2:row+1
    for j = 2:col+1
        temp = A_vert(i-k:i+k, j-k:j+k);
        Sum1 = sum(temp(:));
        Binned(i-1,j-1) = Sum1;
    end
end
end


