import streamlit as st

def local_css(file_name):
    with open(file_name) as f:
        st.markdown('<style>{}</style>'.format(f.read()), unsafe_allow_html=True)

local_css("style.css")
def color_and_font(str,color):
    return "<div> <span class='highlight "+color+"'><span class='bold'>"+str+"</span></span></div>"

def lm(pept_length, alpha):
        """
        Compute weights using linear variation model
        :param pept_length: int
                size of the window
        :param alpha: float between 0 (exclusive) and 1 (inclusive)
                edge weight

        :return: list
                a list of weights.
        """
        weight_lst = []
        if pept_length % 2 != 0:
            for idx in range(0, pept_length):
                if idx <= pept_length // 2:
                    weight = alpha + (1 - alpha) * idx / (pept_length // 2)
                    weight = round(weight, 2)
                    weight_lst.append(weight)
                else:
                    weight = 1 - (1 - alpha) * (idx - pept_length // 2) / (pept_length // 2)
                    weight = round(weight, 2)
                    weight_lst.append(weight)
        else:
            for idx in range(0, pept_length):
                if idx < pept_length / 2:
                    weight = alpha + (1 - alpha) * idx / (pept_length / 2 - 1)
                    weight = round(weight, 2)
                    weight_lst.append(weight)
                else:
                    weight = 1 - (1 - alpha) * (idx - pept_length / 2) / (pept_length / 2 - 1)
                    weight = round(weight, 2)
                    weight_lst.append(weight)

        return weight_lst


def calc_hopp_parker(seq, pep_length, alpha, parker_hydrophilic):
    """
    Calculate the hopp-woods score for each peptide using the linear variation model
    :param seq:str
            protein seq in one-letter code
    :param pep_length:int
            size of the window (length of the peptide)
    :param alpha: float
            edge weight, between 0 (exclusive) and 1 (inclusive)
    :return:tuple
            a tuple (averaged hydrophilicity score, peptide seq)
    """

    # Caculate un-corrected score
    aa_lst = list(seq)
    resi_hopp_lst = [parker_hydrophilic[x] for x in aa_lst]

    # Caculate weights
    weight_lst = lm(pep_length, alpha)
    # st.write("Weights used: ", end="")
    # st.write(weight_lst)

    # a dictionary of {peptide_seq: averaged_hopp_score}
    pept_score_dict = {}

    # Calculate corrected score
    for i in range(0, len(resi_hopp_lst) - pep_length + 1):

        pept_score_lst = resi_hopp_lst[i:i + pep_length]
        weighted_pept_score_lst = []

        for score, weight in zip(pept_score_lst, weight_lst):
            weighted_score = score * weight
            weighted_pept_score_lst.append(weighted_score)

        pept_score = sum(weighted_pept_score_lst) / (sum(weight_lst))  # sum of scores averaged over sum of weights
        pept_seq = "".join(aa_lst[i:i + pep_length])
        pept_score_dict[pept_seq] = pept_score

    # key:value pair was switched in the turple to allow sorting by hopp score
    return [(v, k) for k, v in pept_score_dict.items()]


def lm(pept_length, alpha):
        """
        Compute weights using linear variation model
        :param pept_length: int
                size of the window
        :param alpha: float between 0 (exclusive) and 1 (inclusive)
                edge weight

        :return: list
                a list of weights.
        """
        weight_lst = []
        if pept_length % 2 != 0:
            for idx in range(0, pept_length):
                if idx <= pept_length // 2:
                    weight = alpha + (1 - alpha) * idx / (pept_length // 2)
                    weight = round(weight, 2)
                    weight_lst.append(weight)
                else:
                    weight = 1 - (1 - alpha) * (idx - pept_length // 2) / (pept_length // 2)
                    weight = round(weight, 2)
                    weight_lst.append(weight)
        else:
            for idx in range(0, pept_length):
                if idx < pept_length / 2:
                    weight = alpha + (1 - alpha) * idx / (pept_length / 2 - 1)
                    weight = round(weight, 2)
                    weight_lst.append(weight)
                else:
                    weight = 1 - (1 - alpha) * (idx - pept_length / 2) / (pept_length / 2 - 1)
                    weight = round(weight, 2)
                    weight_lst.append(weight)

        return weight_lst


def calc_hopp(seq, pep_length, alpha ,hopp_scores):
    """
    Calculate the hopp-woods score for each peptide using the linear variation model
    :param seq:str
            protein seq in one-letter code
    :param pep_length:int
            size of the window (length of the peptide)
    :param alpha: float
            edge weight, between 0 (exclusive) and 1 (inclusive)
    :return:tuple
            a tuple (averaged hydrophilicity score, peptide seq)
    """

    # Caculate un-corrected score
    aa_lst = list(seq)
    resi_hopp_lst = [hopp_scores[x] for x in aa_lst]

    # Caculate weights
    weight_lst = lm(pep_length, alpha)
    # print("Weights used: ", end="")
    # print(weight_lst)
    # st.write("Weights used: ", end="")
    # st.write(weight_lst)

    # a dictionary of {peptide_seq: averaged_hopp_score}
    pept_score_dict = {}

    # Calculate corrected score
    for i in range(0, len(resi_hopp_lst) - pep_length + 1):

        pept_score_lst = resi_hopp_lst[i:i + pep_length]
        weighted_pept_score_lst = []

        for score, weight in zip(pept_score_lst, weight_lst):
            weighted_score = score * weight
            weighted_pept_score_lst.append(weighted_score)

        pept_score = sum(weighted_pept_score_lst) / (sum(weight_lst))  # sum of scores averaged over sum of weights
        pept_seq = "".join(aa_lst[i:i + pep_length])
        pept_score_dict[pept_seq] = pept_score

    # key:value pair was switched in the turple to allow sorting by hopp score
    return [(v, k) for k, v in pept_score_dict.items()]


def lm_hopp(pept_length, alpha):
    """
    Compute weights using linear variation model
    :param pept_length: int
            size of the window
    :param alpha: float between 0 (exclusive) and 1 (inclusive)
            edge weight

    :return: list
            a list of weights.
    """
    weight_lst = []
    if pept_length % 2 != 0:
        for idx in range(0, pept_length):
            if idx <= pept_length // 2:
                weight = alpha + (1 - alpha) * idx / (pept_length // 2)
                weight = round(weight, 2)
                weight_lst.append(weight)
            else:
                weight = 1 - (1 - alpha) * (idx - pept_length // 2) / (pept_length // 2)
                weight = round(weight, 2)
                weight_lst.append(weight)
    else:
        for idx in range(0, pept_length):
            if idx < pept_length / 2:
                weight = alpha + (1 - alpha) * idx / (pept_length / 2 - 1)
                weight = round(weight, 2)
                weight_lst.append(weight)
            else:
                weight = 1 - (1 - alpha) * (idx - pept_length / 2) / (pept_length / 2 - 1)
                weight = round(weight, 2)
                weight_lst.append(weight)

    return weight_lst


def calc_hopps(seq, pep_length, alpha, hopp_scores):
    """
    Calculate the hopp-woods score for each peptide using the linear variation model
    :param seq:str
            protein seq in one-letter code
    :param pep_length:int
            size of the window (length of the peptide)
    :param alpha: float
            edge weight, between 0 (exclusive) and 1 (inclusive)
    :return:tuple
            a tuple (averaged hydrophilicity score, peptide seq)
    """

    # Caculate un-corrected score
    aa_lst = list(seq)
    resi_hopp_lst = [hopp_scores[x] for x in aa_lst]

    # Caculate weights
    weight_lst = lm_hopp(pep_length, alpha)
    # st.write("Weights used: ", end="")
    # st.write(weight_lst)

    # a dictionary of {peptide_seq: averaged_hopp_score}
    pept_score_dict = {}

    # Calculate corrected score
    for i in range(0, len(resi_hopp_lst) - pep_length + 1):

        pept_score_lst = resi_hopp_lst[i:i + pep_length]
        weighted_pept_score_lst = []

        for score, weight in zip(pept_score_lst, weight_lst):
            weighted_score = score * weight
            weighted_pept_score_lst.append(weighted_score)

        pept_score = sum(weighted_pept_score_lst) / (
            sum(weight_lst))  # sum of scores averaged over sum of weights
        pept_seq = "".join(aa_lst[i:i + pep_length])
        pept_score_dict[pept_seq] = pept_score

    # key:value pair was switched in the turple to allow sorting by hopp score
    return [(v, k) for k, v in pept_score_dict.items()]