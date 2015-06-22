#!/usr/bin/env python
""" A set of data munging functions to getting things done quickly and cleanly """
import numpy as np


def asList(x):
    """ Convert a variable to a list.

    :param x: an object that you would like to convert to a list. Will raise an
        exception if it is iterable but not already a tuple or list.

    """

    if isinstance(x, list):
        return x
    elif isinstance(x, tuple):
        return list(x)
    elif not hasattr(x, '__iter__'):
        return [x, ]
    else:
        raise ValueError('Could not convert {} to a list'.format(str(x)))


def mergeAndDrop(df, flags, left_on, right_on, flagName, keep_in_list=None, keep_logic=None, drop_in_list=None, drop_logic=None, how='all'):
    """ Merge on a set of flags and drop depending on a criteria.

    Arguments:
        :type df: pandas.DataFrame
        :param df: A data frame that you want to merge flags on to.

        :type flags: pandas.DataFrame
        :param flags: A data frame of flags that you want to merge.

        :param str left_on: Column name that you want to merge in your main data frame.

        :param str right_on: Column name that you want to merge in your flags data frame.

        :param list flagName: Column names with flags.

        :param list keep_in_list: List of what flag values to keep.

        :param str keep_logic: A string with the logical opperator for which
            flags to keep. For example, ('== 0', '<= 1', '>= 10', '!= 2').

        :param list drop_in_list: List of what flag values to drop.

        :param str drop_logic: A string with the logical opperator for which
            flags to drop. For example, ('== 0', '<= 1', '>= 10', '!= 2').

        :type how: str|int
        :param how: If 'all' then a row will be kept if all flags return
            True (keep) or False (drop). If 'any', a row will be kept if any of the
            flags were True (keep) or False (drop). If number between 1 and
            100, then this number will be used as a proportional cutoff. For
            example, if how was 50, then a row will be kept or dropped if a
            flag was True (keep) or False (drop) in 50% of the flags.

    Returns:
        :rtype: pd.DataFrame
        :return: A data frame with the rows that were flagged (keep) or without
            the flagged rows (drop).

    """

    # Merge datasets
    ## Re-index first to make merging easier.
    dfI = df.reset_index()
    flagsI = flags.reset_index()
    merged = dfI.merge(flagsI, how='left', left_on=left_on, right_on=right_on)

    # Make sure flagName is a list
    nameList = asList(flagName)

    # Make sure at least 1 and only 1 keep|drop condition is set
    if len(np.array([keep_in_list, keep_logic, drop_in_list, drop_logic]).nonzero()[0]) > 1:
        print "keep_in_list, keep_logic, drop_in_list, drop_logic are mutually exclusive. Please only provide one."
        raise Exception
    elif len(np.array([keep_in_list, keep_logic, drop_in_list, drop_logic]).nonzero()[0]) == 0:
        print "Please provide at least on of keep_in_list, keep_logic, drop_in_list, drop_logic."
        raise Exception

    # Create Boolean Mask using keep|drop conditions
    if keep_in_list:
        print "Keeping flags in given list."

        # Make sure keep_in_list is a list
        flagList = asList(keep_in_list)

        # Create mask where drop_when values are False
        cleanMask = merged[nameList].isin(flagList)

    elif drop_in_list:
        print "Dropping flags in given list."

        # Make sure keep_in_list is a list
        flagList = asList(drop_in_list)

        # Create mask where drop_when values are False
        cleanMask = ~merged[nameList].isin(flagList)

    elif keep_logic:
        print "Keeping flags with {}.".format(keep_logic)
        cleanMask = eval('merged[nameList] ' + keep_logic)

    elif drop_logic:
        print "Dropping flags with {}.".format(drop_logic)
        cleanMask = ~eval('merged[nameList] ' + drop_logic)

    # Return df with specific rows
    if how == 'all':
        return df[cleanMask.values.all(axis=1)]
    elif how == 'any':
        return df[cleanMask.values.any(axis=1)]
    elif how >= 1 and how <= 100:
        rowMargin = cleanMask.sum(axis=1)
        rowTotal = cleanMask.count(axis=1)
        rowProp = rowMargin / rowTotal * 100
        propMask = rowProp >= how
        return df[propMask.values]
    else:
        print '"how" must have a value of "all", "any", or an integer between 1 and 100'
        raise ValueError


def orderDf(df, varList):
    """ Re-order the columns of a dataframe.

    Arguments:
        :type df: pandas.DataFrame
        :param df: A pandas dataframe

        :param list varList: List of column names you want to be placed at the
            front of your data frame.

    Returns:
        :rtype: pandas.DataFrame
        :returns: A pandas DataFrame with re-ordered columns.

    """
    # Create a list of the other columns in df, not including the columns that
    # are being moved.
    otherCols = [x for x in df.columns if x not in varList]

    # Create new data frame with the columns name in varList at the front.
    dfOrder = df[varList + otherCols].copy()

    return dfOrder

if __name__ == '__main__':
    pass
