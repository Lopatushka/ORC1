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
    
def create_mutation_dict(df):
    """
    Creates a list of dictionaries representing protein mutations with their positions and labels.

    Args:
        df (pd.DataFrame): DataFrame containing a column "Protein change" with mutation information.

    Returns:
        List[dict]: A list of dictionaries, each containing:
            - "position" (int): The numeric position of the mutation in the protein sequence.
            - "label" (str): The full mutation label (e.g., "R846W").
    """
    # List of prtein mutations. Drop NAs
    mutations = list(df["Protein change"].dropna())
    print(f"Number of variats is", len(mutations))

    # Split each string into separate mutations and flatten the list
    mutations_names = [mutation.strip() for group in mutations for mutation in group.split(',')]

    # Create a list of dictionaries with position and label
    mutations_with_labels = [
        {"position": int(re.findall(r'\d+', mutation)[0]), "label": mutation}
        for mutation in mutations_names
    ]
    return mutations_with_labels
    
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

def percantage_aa(seq, A):
    length = len(seq)
    count = 0
    for aa in seq:
        if aa == A:
            count += 1
    return 100*count/length