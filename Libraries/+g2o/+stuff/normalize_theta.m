function theta = normalize_theta(theta)

%return

while (theta < -pi)
    theta = theta + 2* pi;
end

while (theta > pi)
    theta = theta - 2* pi;
end