% read the noisy data u0 in BMP format
u0=imread('terminal.bmp');
%nb=load('noisybrain.mat')
%u0 = nb.Kn
u0=im2double(u0);
[M N]=size(u0);

% visualize the image u0 in Matlab (rescaled) 
imagesc(u0); axis image; axis off; colormap(gray);
%spec_img_log = log(1 + u0);
%imshow(spec_img_log,[]); 


%---------------------------------------------------------------------------
% initialize u by u0 (not necessarily) or by a constant, or by a random image
%---------------------------------------------------------------------------

lam = 0.5:0.5:3
%lam = 0.5

for l = 1:length(lam),
vars = {'u','Energy', 'Fidelity', 'Q'};
clear(vars{:});
    
u=u0;
[M,N]=size(u);


%----------------------------------------------- 
%           PARAMETERS  
%-----------------------------------------------

% coefficient of the TV norm (needs to be adapted for each image) 
lambda=lam(l); 

%  space discretization 
h=1.; 

%  number of iterations (depends on the image) 
IterMax=100;

eps=0.0005; % for regularization
alp = 0.05; % gd rate 

%-----------------------------------------------------------
%     BEGIN ITERATIONS IN ITER 
%-----------------------------------------------------------

l/length(lam)*100

for Iter=1:IterMax, 

	%Iter

    for i=2:M-1,
      for j=2:N-1,

	ux=(u(i+1,j)-u(i-1,j))/(2*h);
    uy=(u(i,j+1)-u(i,j-1))/(2*h);
	uxx = (u(i+1,j)-2*u(i,j)+u(i-1,j))/(h^2);
    uyy = (u(i,j+1)-2*u(i,j)+u(i,j-1))/(h^2);
    uxy = (u(i+1,j+1)-u(i-1,j+1)-u(i+1,j-1)+u(i-1,j-1))/(4*h^2);
    div = (uxx*uy^2-2*ux*uy*uxy+uyy*ux^2)/(eps*eps+ux*ux+uy*uy)^(3/2);
    
	u(i,j)= u(i,j) + alp*sqrt(eps*eps+ux*ux+uy*uy)*(div-2*lambda*(u(i,j)-u0(i,j))) ;

      end
    end

    %------------ FREE BOUNDARY CONDITIONS IN u -------------------
	for i=2:M-1,
          u(i,1)=u(i,2);
          u(i,N)=u(i,N-1);
        end

	for j=2:N-1,
          u(1,j)=u(2,j);
          u(M,j)=u(M-1,j);
        end

        u(1,1)=u(2,2);
        u(1,N)=u(2,N-1); 
        u(M,1)=u(M-1,2);
        u(M,N)=u(M-1,N-1);
%----------------------------------------------------------------------

%---------------- END ITERATIONS IN Iter ------------------------------

%%% Compute the discrete energy at each iteration
en=0.0;  
    for i=2:M-1,
      for j=2:N-1,
      ux=(u(i+1,j)-u(i,j))/h;
      uy=(u(i,j+1)-u(i,j))/h;
      fid=(u0(i,j)-u(i,j))*(u0(i,j)-u(i,j));
      en=en+alp*sqrt(eps*eps+ux*ux+uy*uy)+lambda*fid;
      end
    end
%%% END computation of energy 

Energy(Iter)=en; 

end 


% visualize the image in Matlab (re-scaled)
%figure 
%imagesc(u); axis image; axis off; colormap(gray);
imwrite(u, "tv2_gd/terminal/tv2gd_"+lambda+".jpg", "Quality", 100) 

% use export to save the image in ps or eps format 

% Plot the Energy versus iterations 
%figure 
plot(Energy);legend('Energy/Iterations');
saveas(gcf, "tv2_gd/terminal/tv2gd_Energy_"+lambda+".jpg") 

end



