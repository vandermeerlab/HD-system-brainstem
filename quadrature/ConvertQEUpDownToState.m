function state_tsd = ConvertQEUpDownToState(up_down_tsd) %better vector names

%create empty state vector
state_vec = zeros(1, length(up_down_tsd.tvec)); 

vec_1 = up_down_tsd.data(1,:);
vec_2 = up_down_tsd.data(2,:);

%create and assign states
s1 = 1;
s2 = 2;
s3 = 3;
s4 = 4;

for i = 1:length(vec_1)
   
    if (vec_1(i) == 1 & vec_2(i) == 1)
        state_vec(i) = s1;
        %disp('state 1')
    end
    
    %state 2
    if (vec_1(i) == 1 & vec_2(i) == 0)
        state_vec(i) = s2;
        %disp('state 2')
    end
    
    %state 3
    if (vec_1(i) == 0 & vec_2(i) == 1)
        state_vec(i) = s3;
        %disp('state 3')
    end
    
    %state 4
    if (vec_1(i) == 0 & vec_2(i) == 0)
        state_vec(i) = s4;
        %disp('state 4')
    end
    
end

up_down_tsd.data = state_vec;
state_tsd = up_down_tsd;