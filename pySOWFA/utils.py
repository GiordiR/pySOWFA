from numpy import delete

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
        vector_components.append(x_components)
        vector_components.append(y_components)
        vector_components.append(z_components)

    tensor_components = []
    for tensor_field in tensor_fields:
        xx_components = tensor_field + 'xx'
        yy_components = tensor_field + 'yy'
        zz_components = tensor_field + 'zz'
        xy_components = tensor_field + 'xy'
        xz_components = tensor_field + 'xz'
        yz_components = tensor_field + 'yz'
        tensor_components.append(xx_components)
        tensor_components.append(yy_components)
        tensor_components.append(zz_components)
        tensor_components.append(xy_components)
        tensor_components.append(xz_components)
        tensor_components.append(yz_components)

    return vector_components, tensor_components

def cleanDoubleData(dataArray, coordinateArray1, coordinateArray2, coordinateArray3):
    """ Remove duplicate data and retain their indexes"""
    index = []
    for pos in range(1, len(dataArray)):
        if dataArray[pos] == dataArray[pos-1]:
            index.append(pos)

    newData = delete(dataArray, index)
    newCoordinate1 = delete(coordinateArray1, index)
    newCoordinate2 = delete(coordinateArray2, index)
    newCoordinate3 = delete(coordinateArray3, index)

    return newData, newCoordinate1, newCoordinate2, newCoordinate3
