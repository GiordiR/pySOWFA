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


def fieldsComponentsOF():
    """
    List of OpenFOAM fields components

    :return: vector_components
    :return: tensor_components
    """
    scalar_fields, vector_fields, tensor_fields = fieldsOF()

    vector_components = []
    for vector_field in vector_fields:
        x_components = vector_field + 'x'
        y_components = vector_field + 'y'
        z_components = vector_field + 'z'
        vector_components.append(x_components, y_components, z_components)

    tensor_components = []
    for tensor_field in tensor_fields:
        xx_components = vector_field + 'xx'
        yy_components = vector_field + 'yy'
        zz_components = vector_field + 'zz'
        xy_components = vector_field + 'xy'
        xz_components = vector_field + 'xz'
        yz_components = vector_field + 'yz'
        tensor_components.append(xx_components, yy_components, zz_components, xy_components,
                                 xz_components, yz_components)

    return vector_components, tensor_components
