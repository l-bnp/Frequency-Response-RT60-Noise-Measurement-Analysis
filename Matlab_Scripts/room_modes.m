
% Parameter Definition
c = 343;    % Speed of Sound
L = 15.08;      % Room Lenght
W = 9.8;      % Room Width
H = 8.33;      % Room Height

n=1;

for p = 0:1
    for q = 0:1
        for r = 0:1
            
            modes(n, 2) = p;
            modes(n, 3) = q;
            modes(n, 4) = r;
            modes(n, 5) = (c/2)*sqrt((p/L)^2 + (q/W)^2 + (r/H)^2);
            
            n = n+1;
            
        end
    end   
end

% Remove null rows
modes = modes(~all(modes == 0, 2),:);

% Sort in ascending order, based on frequency
modes = sortrows(modes, 5);

for m = 1:size(modes, 1)
    
    modes(m, 1) = m;
    
end
