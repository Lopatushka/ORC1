import re

def generate_list(x):
    '''Generate a list of all numbers in the given intervals
    Input: list of tuples
    Output: list'''
    generated_list = []
    for start, end in x:
        generated_list.extend(range(start, end + 1))  # `end + 1` to include the endpoint
    return generated_list

def rename_condition(x):
    '''
    Input: pd.Series
    Output: pd.Series
    '''
    if x == 'Meier-Gorlin syndrome 1':
        return x
    elif 'Meier-Gorlin syndrome 1' in x:
        return "Meier-Gorlin syndrome 1 probably"
    elif x == 'ORC1-related disorder':
        return x
    elif 'ORC1-related disorder' in x:
        return 'ORC1-related disorder probably'
    elif x == "Inborn genetic diseases":
        return x
    elif 'Inborn genetic diseases' in x:
        return 'Inborn genetic diseases probably'
    elif x == "not provided":
        return "not provided"
    elif "not provided" in x:
        return "not provided"
    elif x == "not specified":
        return "not provided"
    else:
        return x
    
def mutation_analysis(sb, li):
    '''
    Input - pd.Series with list of mutations to test
    Output - 2 lists of mutations: first - within interval of interest (li); second - outside of interval of interest
    '''
    sb.dropna(inplace=True) # deleete empty positions in the 
    sb = list(sb) # convert pd.Series to list
    sb = [string.split(", ") for string in sb] # split strings with contains several mutaions
    sb = [item for sublist in sb for item in sublist] # flatten list of lists
    interval_mutations = []
    no_interval_mutations = []
    for mutation in sb:
        positions = list(map(int, re.findall(r'\d+', mutation)))
        for position in positions:
            if position in li:
                interval_mutations.append(mutation)
            else:
                no_interval_mutations.append(mutation)
    print(f"Total number of mutations is: {len(sb)}")
    print(f"Number of mutations within interaval of interest is {len(interval_mutations)}")
    print(f"Number of mutations outside interaval of interest is {len(no_interval_mutations)}")
    return(interval_mutations, no_interval_mutations)