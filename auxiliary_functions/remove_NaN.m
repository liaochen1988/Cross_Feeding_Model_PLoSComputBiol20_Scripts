function var_without_NaN = remove_NaN(var)
var_without_NaN = var(~isnan(var));
end

