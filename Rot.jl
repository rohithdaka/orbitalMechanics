function Rot(θ,axis)
    RM1 = [1 0 0; 0 cos(θ) sin(θ); 0 -sin(θ) cos(θ)]
    RM2 = [cos(θ) 0 -sin(θ); 0 1 0 ; sin(θ) 0 cos(θ)]
    RM3 = [cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]
    if axis==1
        return RM1
    elseif axis == 2
        return RM2
    elseif axis == 3
        return RM3
    end
     
end