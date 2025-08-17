function [x,x_left, x_middle, x_right] = generate_array_positions(Ns, Na, d1, d2)
    % Ns: total number of antenna elements
    % Na: number of center elements (with Na-1 gaps of d2)
    % d1: spacing between elements in the side regions
    % d2: spacing between elements in the center region
    % x : array containing the coordinates of all elements (centered at 0)

    % Ensure symmetric structure
    if mod(Ns - Na, 2) ~= 0
        error('Ns - Na must be even to maintain symmetry.');
    end

    Nl = (Ns - Na) / 2;  % number of elements on the left side
    Nr = Nl;             % number of elements on the right side

    % Indices for the center elements: n = Nl : Nl+Na-1
    middle_idx = 0 : Na-1;
    middle_center = (Na - 1) / 2;
    x_middle = (middle_idx - middle_center) * d2;

    % Left elements: extend left from the first center element
    x_left = x_middle(1) - d1 * (Nl:-1:1);

    % Right elements: extend right from the last center element
    x_right = x_middle(end) + d1 * (1:Nr);

    % Concatenate all positions
    x = [x_left, x_middle, x_right];
end