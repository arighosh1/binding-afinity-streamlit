# loading all the libraries necessary
from utils import *
from constants import autoimmune_conditions
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import train_test_split
from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy
import warnings
from streamlit_option_menu import option_menu
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from scipy.stats import pearsonr
import math
from math import sqrt  
import warnings

warnings.filterwarnings("ignore")



def intorduction():
    
    st.write(f"""
        # Introduction
            """)
    st.title("Binding Affinity and Antigenicity Prediction Tool (BAAPT)")
    st.write("Prediction tools for binding affinity and antigenicity are currently webserver-based and do not follow a single software design pattern or reference architecture. Therefore, developing universal models and architectures that may be used as templates for the building of these prediction tools is crucial. As a result, this study explains the Binding Affinity and Antigenicity Prediction Tool (BAAPT)'s concept and pattern, as well as how it is applied to various reverse vaccinology systems. ")
    st.title("FAQ")
    st.write("Q: How to prepare dataset for Binding Affinity?")
    # def header(url):
    #  st.markdown(f'<p style="background-color:#0066cc;color:#33ff33;font-size:12px;border-radius:2%;">{url}</p>', unsafe_allow_html=True)
    # st.write("🎈🎈🎈Dataset Create 1 and 2: https://github.com/arighosh1/BindingAffinity_and_Antigency")
    # header("https://github.com/arighosh1/BindingAffinity_and_Antigency")
    st.write("Q: How can we get the sample dataset for Binding Affinity, protein sequence for Antigenicity, and vaccine data for Autoimmunity prediction?")
    st.write("🎈Binding Affinity:")
    st.write("🎈🎈🎈Training and Testing Dataset (created by Developer): protein-ligand.csv and protein-ligand-test.csv https://github.com/arighosh1/BindingAffinity_and_Antigency")
    st.write("🎈🎈🎈Training and Testing Dataset (taken from published research paper): pdbbind-2016-train.csv and pdbbind-2016-test.csv Jones, D., Kim, H., Zhang, X., Zemla, A., Stevenson, G., Bennett, W. F. D., Kirshner, D., Wong, S. E., Lightstone, F. C., & Allen, J. E. (2021). Improved Protein-Ligand Binding Affinity Prediction with Structure-Based Deep Fusion Inference. Journal of Chemical Information and Modeling, 61(4), 1583–1592.") 
    st.write("🎈Antigenicity:")
    st.write("🎈🎈🎈Protein Sequence: https://www.ncbi.nlm.nih.gov/protein/YP_009724390.1?report=fasta")
    st.write("🎈🎈🎈AntiRNA Antibodies: https://www.uniprot.org/blast/uniprot/B20211004F248CABF64506F29A91F8037F07B67D10159C09")
    st.write("🎈Autoimmunity:")
    st.write("🎈🎈🎈VAERS Data: https://vaers.hhs.gov/data/datasets.html")
    st.write("🎈🎈🎈VAERS Data Use Guide: https://vaers.hhs.gov/docs/VAERSDataUseGuide_en_September2021.pdf")
        
def binding_afinity():
    
    st.title("Binding Affinity")
    st.write("The strength of the binding interaction between a single biomolecule (e.g., protein or DNA) and its ligand/binding partner is referred to as binding affinity (e.g., drug or inhibitor). The equilibrium dissociation constant (KD), which is used to evaluate, and rank order the strengths of bimolecular interactions, is commonly used to measure and report binding affinity. The lower the KD value, the greater the ligand's affinity for its target. The higher the KD value, the weaker the attraction and binding of the target molecule and ligand.")
    st.write("Non-covalent intermolecular interactions such as hydrogen bonding, electrostatic interactions, hydrophobic and Van der Waals forces between the two molecules all influence binding affinity. Furthermore, the presence of other molecules may affect the binding affinity of a ligand to its target molecule.")
    col1, col2 = st.columns([1,1])
    col1.image("before_method.jpg")
    # col2.write("Method: ")
    col2.image("after_method.jpg")
    
    
    #option_menu bar for selection of particular algorithm
    selected2 = option_menu(None, [ "Random Forest","KNN"], 
    icons=['graph-up-arrow', 'graph-up-arrow'], 
    menu_icon="cast", default_index=0, orientation="horizontal")
    
    if selected2 == "KNN":

        st.write("KNN (K — Nearest Neighbors) is one of many (supervised learning) algorithms used in data mining and machine learning; it is a classifier algorithm in which the learning is based on 'how similar' one data (a vector) is to another.")
    
    if selected2 == "Random Forest":
        
        st.write("The random forest, as the name implies, is made up of many individual decision trees that work together as an ensemble. Each individual tree in the random forest produces a class prediction, and the class with the most votes become the prediction of our model.")
    

    my_1 = color_and_font("Upload File Containing PDB_ID(Training Set) : ","blue")
    
    st.markdown(my_1, unsafe_allow_html=True)
    file1 = st.file_uploader("", accept_multiple_files=False)
   
    st.markdown(my_1, unsafe_allow_html=True)
    file2 = st.file_uploader(" ", accept_multiple_files=False)
    if file1 != None and file2!=None :

        df = pd.read_csv(file1)
        df2 = pd.read_csv(file2)

        # Read the data
        # Traning Sets
        X = df.iloc[:, [1, 2]]
        Y = df2.iloc[:, 1]

        X_train, X_valid, Y_train, Y_valid = train_test_split(X, Y, test_size=0.2, random_state=123456)

        # # Optimized parameters
        # ## max_features = 'auto'
        # ## n_estimators=100
        # ## random_state = 1234

        
        # Random Forest
        models_RF_train = {"RF": RandomForestRegressor(bootstrap=True, criterion='friedman_mse', max_depth=None,
                                                       max_features='auto', max_leaf_nodes=None,
                                                       min_impurity_decrease=0.0,
                                                       min_samples_leaf=1, min_samples_split=2,
                                                       min_weight_fraction_leaf=0.0, n_estimators=100,
                                                       n_jobs=None, oob_score=False, random_state=1234,
                                                       verbose=0, warm_start=False)}


        # Calculate the Training and Validation (Refined set) statistics
        scores = {}
        for m in models_RF_train:
            models_RF_train[m].fit(X_train, Y_train)
            scores[m + "_train_r2"] = models_RF_train[m].score(X_train, Y_train)
            Y_pred_valid_rf = models_RF_train[m].predict(X_valid)
            Y_pred_train_rf = models_RF_train[m].predict(X_train)
            scores[m + "_rmse_train"] = sqrt(mean_squared_error(Y_train, Y_pred_train_rf))
            scores[m + "_mae_train"] = mean_absolute_error(Y_train, Y_pred_train_rf)
            scores[m + "_pcc_train"] = pearsonr(Y_train, Y_pred_train_rf)
            scores[m + "_valid_r2"] = r2_score(Y_valid, Y_pred_valid_rf)
            scores[m + "_rmse_valid"] = sqrt(mean_squared_error(Y_valid, Y_pred_valid_rf))
            scores[m + "_mae_valid"] = mean_absolute_error(Y_valid, Y_pred_valid_rf)
            scores[m + "_pcc_valid"] = pearsonr(Y_valid, Y_pred_valid_rf)
        
        if selected2 == "Random Forest":
            st.title("Random Forest Predicted Values")
            scores_RF_train = pd.DataFrame(scores).T
            st.dataframe(scores_RF_train)

        # Calculate statistics for test set (Core set) based on RF model
        scores = {}
        for m in models_RF_train:
            Y_pred_test_rf = models_RF_train[m].predict(X)
            scores[m + "_test_r2"] = r2_score(Y, Y_pred_test_rf)
            scores[m + "_rmse_test"] = sqrt(mean_squared_error(Y, Y_pred_test_rf))
            scores[m + "_mae_test"] = mean_absolute_error(Y, Y_pred_test_rf)
            scores[m + "_pcc_test"] = pearsonr(Y, Y_pred_test_rf)

        if selected2 == "Random Forest":
            scores_RF_test = pd.DataFrame(scores).T
            st.dataframe(scores_RF_test)



        # Save the test prediction result
        Pred_y = pd.DataFrame({'Y_pred_rf': Y_pred_test_rf})
        Exp_y = pd.DataFrame(Y)
        Prediction = pd.concat([Exp_y, Pred_y], axis=1)

  

        YV_array = np.array(Y_valid)
        YT_array = np.array(Y_train)
        XV_array = np.array(X_valid)
        XT_array = np.array(X_train)

        #KNN

        from sklearn.neighbors import KNeighborsRegressor

        knn_model = KNeighborsRegressor(n_neighbors=3)


        knn_model.fit(XT_array, YT_array)

        train_preds = knn_model.predict(XT_array)

        # In[18]:

        # print("KNN Predicted Values:", train_preds)
        if selected2 == "KNN":
            st.title("KNN Predicted Values:")
            st.write( train_preds)

            mse = mean_squared_error(YT_array, train_preds)
            rmse = sqrt(mse)
                    
            st.title("RMSE_train KNN:",str(rmse))

       

        from matplotlib.colors import ListedColormap
        import seaborn as sns

        if selected2 == "KNN":
            #plt.plot(train_preds)
            st.line_chart(train_preds)

def antigencity():
    

    # Amino Acid Scale Defined by Hopp-Woods's original paper
    parker_hydrophilic = {
        "R": 3,
        "D": 3,
        "E": 3,
        "K": 3,
        "S": 0.3,
        "N": 0.2,
        "Q": 0.2,
        "G": 0,
        "P": 0,
        "T": -0.4,
        "A": -0.5,
        "H": -0.5,
        "C": -1,
        "M": -1.3,
        "V": -1.5,
        "I": -1.8,
        "L": -1.8,
        "Y": -2.3,
        "F": -2.5,
        "W": -3.4
    }

    # ## Examples of Usage
    # ### Example 1: Compute Hopp-Woods Scores Without Weights (window=7, $\alpha=1$)
   
    st.title("Antigenicity")
    st.write(
        "The location of continuous epitopes has been linked to properties of polypeptide chains such as hydrophilicity, flexibility, accessibility, turns, exposed surface, polarity, and antigenic propensity. This has resulted in a search for empirical rules that would allow the position of continuous epitopes to be predicted based on specific features of the protein sequence. The propensity scales for each of the 20 amino acids are used in all prediction calculations. Each scale is made up of 20 values that are assigned to each amino acid residue based on their relative proclivity to possess the property described by the scale.")
    st.write("Method:")
    st.write(
        "The amino acids in an interval of the chosen length, centered around residue I are considered when computing the score for a given residue i. In other words, for a window size n, the score for residue I was computed by subtracting the I - (n-1)/2 neighboring residues on each side of the residue. The score for residue I is the average of the scale values for these amino acids, unless otherwise specified (see table 1 for specific method implementation details). In general, a window size of 5 to 7 is appropriate for detecting potentially antigenic regions.")
    st.write("Hopp-Woods scale:")
    st.write(
        "Hopp-Woods scale: Hopp and Woods created their hydrophobicity scale to help them identify potentially antigenic sites in proteins. This scale is essentially a hydrophilic index with negative values assigned to apolar residues. When a window size of 7 is used, antigenic sites are more likely to be predicted.")
    st.write("Parker scale:")
    st.write(
        "In this method, a hydrophilic scale based on peptide retention times was constructed using high-performance liquid chromatography (HPLC) on a reversed-phase column. The epitope region was examined using a seven-residue window. The corresponding scale value was introduced for each of the seven residues, and the arithmetic mean of the seven residue values was assigned to the segment's fourth, (i+3), residue.")
    st.write("1.	Enter a protein sequence in plain format")
    st.write("2.	Select a prediction method")
    st.write("3.	Press Enter")



    # Amino Acid Scale Defined by Hopp-Woods's original paper
    hopp_scores = {
        "R": 0.87,
        "D": 2.46,
        "E": 1.86,
        "K": 1.28,
        "S": 1.50,
        "N": 1.64,
        "Q": 1.37,
        "G": 1.28,
        "P": 0.30,
        "T": 1.15,
        "A": 0.03,
        "H": 0.30,
        "C": 0.11,
        "M": -1.41,
        "V": -1.27,
        "I": -2.45,
        "L": -2.78,
        "Y": -0.78,
        "F": -2.78,
        "W": -3.00
    }


    
    # provide protein seq
    # protein = input("Enter Protein Sequence ")
    protein = st.text_input("Enter Protein Sequence ")
    
    # creating 

    selected3 = option_menu(None, [ "Parker Algorithm","Hopp Algorithm"], 
    icons=['graph-up-arrow', 'graph-up-arrow'], 
    menu_icon="cast", default_index=0, orientation="horizontal")
    
    protein_hopp=protein
    if len(protein) > 0:
        
        
        # calculate averaged Hopp score
        result = calc_hopp(protein, 7, alpha=1, hopp_scores = hopp_scores)

        avg_score = sorted(result, reverse=True)       

        # Plot desired range to show on the x axis.
        # Recommend to change starting position to 1 instead of 0
        x = range(1, 24)

        # range of averaged hopp scores to show on y axis.
        y = [x[0] for x in result[0:23]]

        # plot chart
        if selected3 == "Parker Algorithm":
            st.title("Parker")
            st.title("(Avg Parker Score Sorted, Peptide)")
            st.dataframe(avg_score)
            
            y_expassy=[-0.086, 0.414, 0.086, -0.300, 0.271, 0.271, -0.014, -0.300,
             -0.800, -0.543, -0.329, -1.014,  -1.057 , -0.943, -0.657,
              -0.843, -0.343, -0.343, -0.043, -0.000, 0.171, 0.086,0.343,
            ]
            
            st.title("Parker Graph")
            plt.figure(1)
            plt.plot(x, y, "r-", linewidth=7, alpha=0.4)
            plt.plot(x, y_expassy, "b--")
            plt.xlabel("Amino Acid Position")
            plt .ylabel("Hydrophilicity Score")
            plt.legend(["Script", "Expasy"], loc="lower right")
            plt.savefig("parker",dpi=100)
            st.image("parker.png")
        
        



        # calculate averaged Hopp score
        
        if selected3 == "Hopp Algorithm":
            # hopp start
            hopp_scores = {
                "R": 3,
                "D": 3,
                "E": 3,
                "K": 3,
                "S": 0.3,
                "N": 0.2,
                "Q": 0.2,
                "G": 0,
                "P": 0,
                "T": -0.4,
                "A": -0.5,
                "H": -0.5,
                "C": -1,
                "M": -1.3,
                "V": -1.5,
                "I": -1.8,
                "L": -1.8,
                "Y": -2.3,
                "F": -2.5,
                "W": -3.4
            }
            
            col1 , col2 = st.columns(2)
            pep_length = col1.radio("Window Size",[7,9,7])
            alpha = col2.radio("Alpha value (/=10)",[1,5,1])
            alpha /= 10
              
            result = calc_hopps(protein_hopp, pep_length, alpha, hopp_scores=hopp_scores)

            result_1 = sorted(result, reverse=True)
           
            # Plot desired range to show on the x axis.
            # Recommend to change starting position to 1 instead of 0
            x = range(1, 24)

            # range of averaged hopp scores to show on y axis.
            y = [x[0] for x in result[0:23]]

        # plot chart
        if selected3 == "Hopp Algorithm":
            st.title("Hopp Score")
            st.title("(Avg Hopp Score Sorted, Peptide)")
            st.dataframe(result_1)
            
            # ## Validating against Expasy Result

            y_expassy2 = [-0.176, 0.059, 0.335, 0.379, 0.276, -0.182, -0.250, -0.018,
                        -0.253, -0.632, -0.968, -0.994, -0.932, -0.909, -0.921, -0.738,
                        -0.618, -0.247, -0.097, 0.344, 0.221, 0.256, 0.115
                        ]
            
            st.title("Hopp Graph")
            plt.figure(2)
            plt.plot(x, y, "r-", linewidth=7, alpha=0.4)
            plt.plot(x, y_expassy2, "b--")
            
            plt.xlabel("Amino Acid Position")
            plt.ylabel("Hydrophilicity Score")
            plt.legend(["Script", "Expasy"], loc="lower right")
            plt.savefig("hopp", dpi=100)
            st.write("Window Size = ", pep_length)
            st.write("Alpha Value = ", alpha)
            st.image("hopp.png")

        

def autoimmunity():
    
    import pandas as pd
    import numpy as np
    from pymatch.Matcher import Matcher
    import re
    from tqdm.autonotebook import tqdm
    import warnings
    from IPython.testing.globalipapp import get_ipython
    import streamlit as st
    st.title("Autoimmunity Prediction for COVID-19")
    st.write("This study looks at the long-term consequences of immunization. It investigates a group of unfavorable occurrences known as 'autoimmune diseases,' which involve the immune system attacking its own cells and include LLD, rheumatoid arthritis, and Crohn's disease. The following types of adverse occurrences are considered. Autoimmunity is a serious issue for COVID-19 vaccination recipients. Autoimmune illnesses are unpleasant, long-lasting, and can reduce a person's earning potential. Certain autoimmune diseases can be deadly or severe in rare circumstances. We used the Vaccine Adversity Reporting System (VAERS) input data set (1990-2021). It's a countrywide early warning system for possible security concerns with American vaccinations. The results demonstrate that the COVID-19 vaccination is safe for those with autoimmune diseases. When a person reports a safety signal, the risk of receiving a COVID-19 vaccination is only twice as high as when a person receives a non-COVID-19 vaccine. It's important to note that the link isn't causative. These investigations do not establish a cause; rather, they aid in the identification of possible issues that may be investigated further using other approaches.")

    tqdm.pandas()


    warnings.filterwarnings('ignore')

    ip=get_ipython()
    ip.run_line_magic('matplotlib', 'inline')

    # Fixed seed for reproducibility
    np.random.seed(4072021)

    # print(f"Python version: {sys.version}")
    # print(f"OS version: {platform.platform()}")
    # print(f"pandas version: {pd.__version__}")
    # print(f"numpy version: {np.__version__}")
    # print(f"scipy version: {scipy.__version__}")
    # print(f"statsmodels version: {statsmodels.__version__}")


    # In[4]:

    vac = st.file_uploader("Enter your (Vaccine information (e.g., vaccine name, manufacturer, lot number, route, site, and number of previous doses administered)) : ")
    rec = st.file_uploader("Enter your (Other vaccine information) : ")
    sym = st.file_uploader("Enter your (Symptoms information) : ")
    st.write("🎈Remember that positive predictive value (PPV) varies with disease prevalence when interpreting results from diagnostic tests. Prevalence is a measure of how common a disease is in a specified at-risk population at a specific time point or period. PPV is the percent of positive test results that are true positives. As disease prevalence decreases, the percent of test results that are false positives increase.")
    st.write("🎈For example, a test with 98% specificity would have a PPV of just over 80% in a population with 10% prevalence, meaning 20 out of 100 positive results would be false positives.")
    st.write("🎈The same test would only have a PPV of approximately 30% in a population with 1% prevalence, meaning 70 out of 100 positive results would be false positives.  This means that, in a population with 1% prevalence, only 30% of individuals with positive test results actually have the disease.")
    st.write("🎈At 0.1% prevalence, the PPV would only be 4%, meaning that 96 out of 100 positive results would be false positives.")
    st.write("🎈Health care providers should take the local prevalence into consideration when interpreting diagnostic test results.")
    st.write("🎈Consider positive results in combination with clinical observations, patient history, and epidemiological information.")
    st.write("🎈A false negative is a test result that is wrong, because it indicates the person is not infected when they really are, or that they don’t have antibodies when they actually do.")
    st.write("🎈A false positive is a test result that is wrong, because it indicates the person is infected when they really are not or that they have antibodies when they really don’t.")
    if(vac!=None and rec!=None and sym!=None):
        vax_frames = []

        df = pd.read_csv(vac, index_col=None, header=0, encoding="latin")
        vax_frames.append(df)

        vax = pd.concat(vax_frames, axis=0, ignore_index=True)[["VAERS_ID", "VAX_TYPE"]]
        vax["VAX_TYPE"] = vax["VAX_TYPE"] == "COVID19"
        vax.columns = ["VAERS_ID", "IS_COVID_VACCINE"]


        # In[3]:




        recipient_frames = []

        df = pd.read_csv(rec, index_col=None, header=0, encoding="latin")
        recipient_frames.append(df)

        recipients = pd.concat(recipient_frames, axis=0, ignore_index=True)[["VAERS_ID", "SEX", "CAGE_YR"]]


        # In[5]:


        age_bands = {0: "<18",
                     18: "18-25",
                     26: "26-40",
                     41: "41-55",
                     56: "56-70",
                     71: ">70",
                     int(recipients.CAGE_YR.max()): "max"}

        recipients["AGE"] = pd.cut(recipients.CAGE_YR, bins=list(age_bands.keys()), labels=list(age_bands.values())[:-1])
        recipients = recipients.drop("CAGE_YR", axis=1).dropna()


        # In[6]:



        symptoms_frames = []

        df = pd.read_csv(sym, index_col=None, header=0)
        symptoms_frames.append(df)


        symptoms = pd.melt(pd.concat(symptoms_frames, axis=0, ignore_index=True)[["VAERS_ID", "SYMPTOM1", "SYMPTOM2", "SYMPTOM3", "SYMPTOM4", "SYMPTOM5"]],
                       id_vars="VAERS_ID",
                       value_vars=(f"SYMPTOM{i}" for i in range(1, 6))).drop("variable", axis=1)

        symptoms.columns = ("VAERS_ID", "SYMPTOM")




        vaccination_data = vax.merge(recipients, how="inner", on="VAERS_ID")


        p_normals = r".*negative$|.*\snormal$|.*(scopy|graphy|gram|metry|opsy)$|.*(count|percentage|level|test|assay|culture|X-ray|imaging|gradient|band(s)?|index|surface area|gas|scale|antibod(y|ies)|urine absent|Carotid pulse|partial pressure|time|P(C)?O2)$|Oxygen saturation$|End-tidal.*"
        p_tests = r".*(ase|ose|ine|enzyme|in|ine|ines|ium|ol|ole|ate|lytes|ogen|gases|oids|ide|one|an|copper|iron)$|.*(level therapeutic)$|.*(globulin)\s.{1,2}$|Barium (swallow|enema)"
        p_procedures = r".*(plasty|insertion|tomy|ery|puncture|therapy|treatment|tripsy|operation|repair|procedure|bypass|insertion|removal|graft|closure|implant|lavage|support|transplant|match|bridement|application|ablation)$|Incisional drainage$|.* stimulation$|Immunisation$"
        p_normal_procedures = r"(Biopsy|pH|.* examination|X-ray|.* pulse|Blood|Electro(.*)gram|.* test(s)?|Echo(.*)gram|.*(scopy)|Cardiac (imaging|monitoring|ventriculogram)|Chromosomal|Carbohydrate antigen|Cell marker|.* examination|Computerised tomogram|Culture|.* evoked potential(s)?|Cytology|Doppler)(?!.*(abnormal|increased|decreased|depression|elevation|present|absent))"
        p_managements = r"(Catheter|Device\).*|.* care$|.* user$|Cardiac pacemaker .*"
        p_other_irrelevants = r"Blood group.*|Blood don(or|ation)$|Drug (abuse(r)?|dependence|screen).*|Elderly|Non-tobacco user|No adverse event"
        p_covid_related = r".*COVID-19(prophylaxis|immunisation|screening)|Asymptomatic COVID-19"

        p = re.compile("|".join([p_normals, p_tests, p_procedures, p_normal_procedures, p_other_irrelevants, p_covid_related]))



        symptoms = symptoms[symptoms.SYMPTOM.str.match(p) == False]


        symptoms["IS_AUTOIMMUNE"] = symptoms.SYMPTOM.isin(autoimmune_conditions)


        instances = vax.merge(symptoms[["VAERS_ID", "IS_AUTOIMMUNE"]].groupby("VAERS_ID").agg({"IS_AUTOIMMUNE": np.max}), how="inner", left_on="VAERS_ID", right_index=True).merge(recipients, how="inner").groupby("VAERS_ID").first()


        maj_min_value=pd.crosstab(instances.IS_COVID_VACCINE, instances.IS_AUTOIMMUNE)
        st.title(maj_min_value)
        st.write(maj_min_value)

        cases = instances[instances.IS_AUTOIMMUNE == True]
        controls = instances[instances.IS_AUTOIMMUNE == False]

        cases["VAERS_ID"] = cases.index
        controls["VAERS_ID"] = controls.index

        m = Matcher(cases, controls, yvar="IS_AUTOIMMUNE", exclude=["VAERS_ID"])