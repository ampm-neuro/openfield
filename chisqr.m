function chisqr(observed, expected)

chi_value = sum(sum(((observed - expected).^2)./expected))

if min(size(observed)) > 1

    degrees_freedom = (length(observed(:,1)) - 1) * (length(observed(1,:)) - 1)
    
else
    
    degrees_freedom = max(size(observed)) - 1
    
end

p_value = 1 - chi2cdf(chi_value,degrees_freedom)

end

