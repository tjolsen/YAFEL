function [  ] = cartesian_FMM()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = 50;
x=1:N;
y=1:N;
[X,Y] = meshgrid(x,y);
frozen = zeros(N);
phi = zeros(N);

%initialize circle of radius 10 at (50, 50)
frozen(1,1) = 1;
phi(1,1) = -10;


% keyboard;
count = 0;
trialList = [];
for outerI = 1:N
    outerI
    for outerJ = 1:N
        if(frozen(outerI, outerJ) ~= 1)
            continue;
        end
        
        trialList = [outerI, outerJ, phi(outerI, outerJ)];
        
        do = 1;
        while(size(trialList,1)>0 || do==1)
            do=0;

%             size(trialList)
            r = getMin(trialList,3);
            ii = trialList(r,1);
            jj = trialList(r,2);
            key = trialList(r,3);
            
            trialList = trialList(1:size(trialList,1) ~= r, :);
            
            frozen(ii,jj) = 1;
            phi(ii,jj) = key;
            
            
            x = X(ii,jj);
            y = Y(ii,jj);
            
            % Get normal direction of front
            grad = [0 0];
            singular = 0;
            singularx = 0;
            singulary = 0;
            if(ii>1 && frozen(ii-1,jj))
                grad(2) = (phi(ii,jj) - phi(ii-1,jj))/(Y(ii,jj)-Y(ii,jj));
            elseif(ii<N && frozen(ii+1,jj))
                grad(2) = (phi(ii+1,jj)-phi(ii,jj))/(Y(ii+1,jj)-Y(ii,jj));
            else
                singulary = 1;
            end
            
            if(jj>1 && frozen(ii,jj-1))
                grad(1) = (phi(ii,jj) - phi(ii,jj-1))/(X(ii,jj)-X(ii,jj-1));
            elseif(jj<N && frozen(ii,jj+1))
                grad(1) = (phi(ii,jj+1) - phi(ii,jj))/(X(ii,jj+1)-X(ii,jj));
            end
            
            %Check if gradient of phi is singular in one or both directions
            TOL = 1e-6;
            if(norm(grad) < TOL)
                singular = 1;
            end
            
            if(abs(grad(1)) < TOL)
                singularx = 1;
            end
            
            if(abs(grad(2)) < TOL)
                singulary = 1;
            end
            
            n = grad/norm(grad);
            
            % Add neighbors to trialList
            newEntries = [];
            if(ii > 1)
                if(~frozen(ii-1,jj))
                    dy = y-Y(ii-1,jj);
                    if(singulary)
                        n(2) = 1;
                    end
                    newkey = key + n(2)*dy;
                    newEntries = [newEntries; ii-1, jj, newkey];%, x, y, key];
                end
            end
            if(ii < N)
                if(~frozen(ii+1,jj))
                    dy = Y(ii+1,jj)-y;
                    if(singulary)
                        n(2) = 1;
                    end
                    newkey = key + n(2)*dy;
                    newEntries = [newEntries; ii+1, jj, newkey];%, x, y, key];
                end
            end
            if(jj > 1)
                if(~frozen(ii,jj-1))
                    dx = x - X(ii,jj-1);
                    if(singularx)
                        n(1) = 1;
                    end
                    newkey = key + n(1)*dx;
                    newEntries = [newEntries; ii, jj-1, newkey];%, x, y, key];
                end
            end
            if(jj < N)
                if(~frozen(ii,jj+1))
                    dx = X(ii,jj+1)-x;
                    if(singularx)
                        n(1) = 1;
                    end
                    newkey = key + n(1)*dx;
                    newEntries = [newEntries; ii, jj+1, newkey];%, x, y, key];
                end 
            end
            
            trialList = [trialList; newEntries];
            contourf(X,Y,phi)
        end
        
%         contourf(X,Y,phi)
%         drawnow
    end
end

keyboard
end

function r = getMin(list,col)
TOL= 1e-5;

[r,c] = find(abs(list(:,col)-min(list(:,col)))<TOL,1,'first');
end