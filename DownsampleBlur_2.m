%% This function takes a datacube A as input and blurs each of the frames by 
%  multiplying it with the blurring matrix PSF

function C=DownsampleBlur_2(A,PSF)
nim=size(A,3);
dx = size(A,1);
dy = size(A,2);
for t=1:nim
    temp = (vec(A(:,:,t)));
    Blur_vec = PSF*temp;
    Blur = reshape(Blur_vec, [dx,dy]);
    Blur_dwnsmp = Blur(1:4:dx,1:4:dy);
    C(:,:,t) = Blur_dwnsmp; 
end

end