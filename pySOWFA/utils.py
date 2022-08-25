def fieldsOF():
    """
    List of OpenFOAM fields

    :return: vector_fields
    :return: scalar_fields
    :return: tensor_fields
    """
    scalar_fields = ['p', 'pMean', 'nuSgs', 'k', 'kMean', 'ppMean']
    vector_fields = ['U', 'UMean', 'omega', 'upMean']
    tensor_fields = ['RMean', 'uuMean', 'uuRTotal', 'UPrime2Mean']

    return scalar_fields, vector_fields, tensor_fields
