import streamlit as st
from streamlit_option_menu import option_menu
import pagesnew as pg
from PIL import Image

img = Image.open('icon1.jpg')
st.set_page_config(layout="wide",
                page_title='RATIONAL VACCINE DESIGN FOR VIRUS USING MACHINE LEARNING APPROACHES', 
                page_icon=img,
                initial_sidebar_state='auto')



st.markdown("""
<style>
.app-header {
    font-size:50px;
    color: #F63366;
    font-weight: 700;
}
.sidebar-header{
    font-family: "Lucida Sans Unicode", "Lucida Grande", sans-serif;
    font-size: 28px;
    letter-spacing: -1.2px;
    word-spacing: 2px;
    color: #FFFFFF;
    font-weight: 700;
    text-decoration: none;
    font-style: normal;
    font-variant: normal;
    text-transform: capitalize;
}
.positive {
    color: #000000;
    font-size:30px;
    font-weight: 700;  
}
.negative {
    color: #70F140;
    font-size:30px;
    font-weight: 700;  
}
</style>
""", unsafe_allow_html=True)


st.markdown(
    '''
        <style>
            @media (max-width: 991.98px)
            {
                .sidebar .sidebar-content 
                {
                    background-color: #003366;
                }
            }
        </style>
    ''',
    unsafe_allow_html=True
)



with st.sidebar:
    col1,col2 = st.columns([2,4])
    col1.image("https://github.com/arighosh1/BindingAffinity_and_Antigency/blob/main/icon1.jpg?raw=true",width=300,output_format='PNG')
    
    selected = option_menu(
        menu_title = "",
        icons = ["grid-fill","calendar","graph-up-arrow","book",'gear'],
        menu_icon = 'fire',
        options = ["Introduction","Binding Affinity",'Antigencity',"Autoimmunity","Setting"],
        default_index = 0
    )
        
              
if selected == "Introduction":
    pg.intorduction()
if selected == "Binding Affinity":
    pg.binding_afinity()
if selected =="Antigencity":
    pg.antigencity()
if selected =="Autoimmunity":
    pg.autoimmunity()
        