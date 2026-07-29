#MUnkle Molecular Unknown Puzzle
#© 2026 Kelsey Sakimoto

import streamlit as st #to run streamlit webapp
from streamlit_ketcher import st_ketcher #for structure editor
from streamlit_local_storage import LocalStorage #for user play history
from pandas import read_csv #for csv import
from plotly.subplots import make_subplots #for plotting spectra
import plotly.graph_objects as go

#for chemoinfomatics and structure drawing
from rdkit import Chem
from rdkit.Chem import Draw,rdRascalMCES
from rdkit.Chem.Draw import rdMolDraw2D

st.iframe(("https://organicchemistrydata.org/js/ketcher/index.html?hiddenControls="
           "arom,dearom,cip,check,analyse,recognize,miew,reaction-plus,arrows,"
           "reaction-mapping-tools,rgroup,shape,text,images,enhanced-stereo,fullscreen,"
           "sgroup,open,save,paste,copy,cut,copy-mol,create-monomer,shape-ellipse,"
           "shape-rectangle,shape-line,extended-table,any-atom,structure-library"))

#global display options
structure_width = 300
structure_height = 200

#create local storage object
local_storage = LocalStorage()

#setup play history file if not already present in local storage
if 'history' not in local_storage.getAll():
    local_storage.setItem('history', ['','','',''])

#create popup to clear play history
@st.dialog('Are you sure?')
def clearRecord():
    st.write('This will erase your wins and losses for all MUnkles.')
    if st.button('Clear play history'):
        local_storage.deleteAll()
        st.rerun()

#create session variable on first play or refresh
if 'puzzNum' not in st.session_state:
    st.session_state['puzzNum'] = 1 #tracks current puzzle
    st.session_state.lastPuzz = 0 #previous puzzle, to cut down on reruns
    st.session_state.last_highC = False #previous high contrast state
    st.session_state.last_expH = False #previous explicit H state
    st.session_state.ketcher_version = 0 #to rerender ketcher on new puzzles

#function to load answer key datafiles, create plots as cached function to speed up app
@st.cache_data(show_spinner='Loading MUnkle')
def setupKey(puzz):
    key_file = str(puzz) + '/' + str(puzz) + '.txt' #molecule information file

    keyInfo = read_csv(key_file,sep = ':',index_col=0,skipinitialspace=True) #dict will contain 'Name', 'SMILES', 'Info', 'HNMR' and 'CNMR'
    canonKeySMILES = Chem.CanonSmiles(keyInfo.at['SMILES','Value']) #canonize SMILES string
    molKey = Chem.MolFromSmiles(canonKeySMILES) #get Mol object of answer 
    molKeyH =Chem.AddHs(Chem.MolFromSmiles(canonKeySMILES)) #get version of mol object with explicit Hs

    spec_files = {'HNMR' : str(puzz) + '/' + str(puzz) + '_HNMR.csv',
                  'CNMR' : str(puzz) + '/' + str(puzz) + '_CNMR.csv',
                  'IR' : str(puzz) + '/' + str(puzz) + '_IR.csv',
                  'MS' : str(puzz) + '/' + str(puzz) + '_MS.csv'}
    spectra = {}
    
    for spec in spec_files:
        data = read_csv(spec_files[spec])
        spectra[spec] = data

    d = rdMolDraw2D.MolDraw2DSVG(structure_width,structure_height) #to draw molecules
    dopts = d.drawOptions()
    dopts.useBWAtomPalette()

    dH = rdMolDraw2D.MolDraw2DSVG(structure_width,structure_height) #to draw molecules
    dHopts = dH.drawOptions()
    dHopts.useBWAtomPalette()
    
    d.DrawMolecule(molKey)
    d.FinishDrawing()
    keyImage = d.GetDrawingText()

    dH.DrawMolecule(molKeyH)
    dH.FinishDrawing()
    keyImageH = dH.GetDrawingText()

    return keyInfo, canonKeySMILES, molKey, molKeyH, spectra, keyImage,keyImageH

#to setup a new puzzle
def newPuzzle():
    if st.session_state.puzzNum != st.session_state.lastPuzz: #checks if it's a new puzzle on rerun
        st.session_state.guessCounter = 0 #reset guess counter
        puzz = st.session_state['puzzNum']
        
        st.session_state.key, st.session_state.canonKeySMILES, st.session_state.molKey, st.session_state.molKeyH, st.session_state.spectra, st.session_state.imageKey, st.session_state.imageKeyH = setupKey(puzz)

        st.session_state.guessSMILES = []
        st.session_state.guessImage = []
        st.session_state.win = False
        st.session_state.lastPuzz = puzz  

@st.cache_data(show_spinner = 'Checking guess...')
def checkGuess(molKey,molGuess,highContrast):
    opts = rdRascalMCES.RascalOptions()
    opts.allBestMCESs = True
    opts.similarityThreshold = 0

    results = rdRascalMCES.FindMCES(molKey, molGuess, opts)

    mostSimilar = None
    mostSimilarSize = 0

    for res in results:
        if res.similarity > mostSimilarSize:
            mostSimilarSize = res.similarity
            mostSimilar = res
    
    bondHighlights_max = {}
    bondHighlights_rest = {}
    atomHighlights_max = {}
    atomHighlights_rest = {}
    remainderHighlights = {}
    
    #get highlights for all matching atoms/bonds, color as "correct, but not in the right place"
    matchedAtomsKey = []
    matchedAtomsGuess = []
    remainderKey = []
    remainderGuess = []
    remainderMatches = []
    
    correctFragColor = (200/255, 182/255, 83/255)
    if highContrast:
        correctFragColor = (130/255, 170/255, 215/255)

    maxFragColor = (108/255, 169/255, 101/255)
    if highContrast:
        maxFragColor = (245/255, 120/255, 90/255)

    if mostSimilar != None:
        for bondPair in mostSimilar.bondMatches():
            bondHighlights_rest[bondPair[1]] = correctFragColor
        
        for atomPair in mostSimilar.atomMatches():
            atomHighlights_rest[atomPair[1]] = correctFragColor
            matchedAtomsKey.extend([atomPair[0]])
            matchedAtomsGuess.extend([atomPair[1]])
  
    for atom in range(molKey.GetNumAtoms()):
        if atom not in matchedAtomsKey:
            remainderKey += [atom]
    
    for atom in range(molGuess.GetNumAtoms()):
        if atom not in matchedAtomsGuess:
            remainderGuess += [atom]
    
    for atomIndexKey in remainderKey:
        for atomIndexGuess in remainderGuess:
            if molKey.GetAtomWithIdx(atomIndexKey).GetSmarts() == molGuess.GetAtomWithIdx(atomIndexGuess).GetSmarts():
                remainderMatches += [(atomIndexKey,atomIndexGuess)]
                remainderGuess.remove(atomIndexGuess)
   
    for atomPair in remainderMatches:
        remainderHighlights[atomPair[1]] = correctFragColor

    if mostSimilar != None:
        mostSimilar.largestFragmentOnly()
        
        for bondPair in mostSimilar.bondMatches():
            bondHighlights_max[bondPair[1]] = maxFragColor
        
        for atomPair in mostSimilar.atomMatches():
            atomHighlights_max[atomPair[1]] = maxFragColor
    
    bondHighlights = bondHighlights_rest | bondHighlights_max
    atomHighlights = atomHighlights_rest | atomHighlights_max | remainderHighlights
    
    return bondHighlights, atomHighlights

st.title(":orange[MUnkle]", text_alignment = "center")
st.header(":orange[M]olecular :orange[Unk]nown Puzz:orange[le]", text_alignment = "center") 
        
headerLeft, headerRight = st.columns(2,vertical_alignment='bottom')

with headerLeft:
    record = local_storage.getItem('history')
    st.selectbox("**Pick a puzzle**",range(1,5), key='puzzNum', on_change=newPuzzle(), width=100,format_func=lambda x: str(x)+ " " + record[x-1])

with headerRight:
    with st.container(horizontal=True,horizontal_alignment='right',vertical_alignment='bottom'):
        with st.popover("",icon=":material/help:"):
            st.text('rules and instructions')
            if st.button('Clear play history?'):
                clearRecord()
        with st.popover("",icon=":material/info:"):
            st.text('stuff about game, credits for spectra')
        with st.popover("",icon=":material/feedback:"):
            st.text('comments and suggestions')
        with st.popover("",icon=":material/volunteer_activism:"):
            st.text('user submissions')
            st.text('link to funding')

chart1, chart2 = st.columns(2)
chart3, chart4 = st.columns(2)

tab1, tab2, tab3, tab4 = st.tabs(['MS','FTIR','¹H NMR','¹³C NMR'],)

with tab1:
    figMS = go.Figure()
    figMS.add_trace(go.Bar(x=st.session_state.spectra['MS']['m/z'], y=st.session_state.spectra['MS']['Intensity'],name='',hovertemplate='%{x} m/z, %{y} rel. int.'))
    
    figMS.update_layout(title=dict(text="MS"), showlegend=False, xaxis=dict(title='m/z'),yaxis=dict(title='Relative Intensity'))
    st.plotly_chart(figMS,width='stretch')
    
with tab2:
    figIR = go.Figure()
    figIR.add_trace(go.Scatter(x=st.session_state.spectra['IR']['Wavenumber (cm-1)'], y=st.session_state.spectra['IR']['Transmittance (%)'],mode='lines',name='ν̃',hovertemplate='%{x} cm⁻¹'))
    
    figIR.update_layout(title=dict(text="FTIR"), showlegend=False, xaxis=dict(autorange="reversed",title='Wavenumber (cm⁻¹)'),yaxis=dict(title='Transmittance (%)',showgrid=False,showticklabels=False))
    st.plotly_chart(figIR,width='stretch')
    
with tab3:
    figHNMR = make_subplots(specs=[[{"secondary_y": True}]])
    
    # Add traces
    figHNMR.add_trace(
        go.Scatter(x=st.session_state.spectra['HNMR']['Chemical Shift (ppm)'], y=st.session_state.spectra['HNMR']['Intensity'],mode='lines',name='δ',hovertemplate='%{x} ppm'),
        secondary_y=True,
    )
    
    figHNMR.add_trace(
        go.Scatter(x=st.session_state.spectra['HNMR']['Chemical Shift (ppm)'], y=st.session_state.spectra['HNMR']['Integral'],mode='lines',name='∫',hovertemplate='%{y}'),
        secondary_y=False,
    )
    figHNMR.update_layout(title=dict(text="¹H NMR", subtitle=dict(text=st.session_state.key.at['HNMR','Value'])), showlegend=False, xaxis=dict(autorange="reversed",title='Chemical Shift (ppm)'), yaxis2=dict(showgrid=False,showticklabels=False),yaxis=dict(title='Integral', ticksuffix=' H'))
    st.plotly_chart(figHNMR,width='stretch')
    
with tab4:
    figCNMR = go.Figure()
    figCNMR.add_trace(go.Scatter(x=st.session_state.spectra['CNMR']['Chemical Shift (ppm)'], y=st.session_state.spectra['CNMR']['Intensity'],mode='lines',name='δ',hovertemplate='%{x} ppm'))
    
    figCNMR.update_layout(title=dict(text="¹³C NMR", subtitle=dict(text=st.session_state.key.at['CNMR','Value'])), showlegend=False, xaxis=dict(autorange="reversed",title='Chemical Shift (ppm)'),yaxis=dict(showgrid=False,showticklabels=False))
    st.plotly_chart(figCNMR,width='stretch')

with st.expander('**Draw a structure to guess a molecule**'):
    ketcher_key = f"ketcher_{st.session_state.puzzNum}"
    currentSMILES = st_ketcher(key=ketcher_key)

col1, col2 = st.columns(2)
with col1:
    highC = st.checkbox("High Contrast Mode?")
with col2:
    expH = st.checkbox("Show explicit hydrogens?")

if Chem.MolFromSmiles(currentSMILES) is None:
    st.subheader("Not a valid structure, try again", text_alignment='center')
    currentSMILES = ""

canonSMILES = Chem.CanonSmiles(currentSMILES)

if canonSMILES != "" and st.session_state.win is False and st.session_state.guessCounter < 6 and highC == st.session_state.last_highC and expH == st.session_state.last_expH:
    if canonSMILES in st.session_state.guessSMILES:
        st.subheader("Structure already guessed", text_alignment='center')

    else:
        st.session_state.guessSMILES.append(canonSMILES)
        molGuess = Chem.MolFromSmiles(canonSMILES)
        molGuessH = Chem.AddHs(Chem.MolFromSmiles(canonSMILES))
        
        highlights = {'noH':{'noHC': checkGuess(st.session_state.molKey,molGuess,False),
                             'HC':checkGuess(st.session_state.molKey,molGuess,True)},
                      'H': {'noHC': checkGuess(st.session_state.molKeyH,molGuessH,False),
                            'HC': checkGuess(st.session_state.molKeyH,molGuessH,True)}}
        image_strings = {}
        
        for h in highlights:
            if h == 'H':
                mol = molGuessH
            else:
                mol = molGuess
                
            for hc in highlights[h]:
                d = rdMolDraw2D.MolDraw2DSVG(structure_width,structure_height) #to draw molecules
                dopts = d.drawOptions()
                dopts.useBWAtomPalette()
                
                d.DrawMolecule(mol, highlightAtoms=highlights[h][hc][1].keys(),highlightBonds=highlights[h][hc][0].keys(),
                               highlightAtomColors=highlights[h][hc][1],highlightBondColors=highlights[h][hc][0])
                d.FinishDrawing()
                image_strings[hc+h] = d.GetDrawingText()

        st.session_state.guessImage.append(image_strings)
        st.session_state.guessCounter +=1

        if canonSMILES == st.session_state.canonKeySMILES:
            st.session_state.win = True

row3 = st.columns(3, border=True)
row4 = st.columns(3, border=True)

for index,cont in enumerate(row3+row4):
    if index > st.session_state.guessCounter-1:
        cont.empty()
    elif highC and expH:
        cont.image(st.session_state.guessImage[index]['HCH']) 
    elif highC and not expH:
        cont.image(st.session_state.guessImage[index]['HCnoH'])
    elif not highC and expH:
        cont.image(st.session_state.guessImage[index]['noHCH']) 
    elif not highC and not expH:
        cont.image(st.session_state.guessImage[index]['noHCnoH'])

    st.session_state.last_highC = highC
    st.session_state.last_expH = expH

if st.session_state.guessCounter >=6 or st.session_state.win:
    with st.container(border = True, horizontal_alignment='center'):
        if st.session_state.win:
            st.header("Molecule Identified!",text_alignment='center')
            st.balloons()
            hist1 = local_storage.getItem('history')
            hist1[st.session_state['puzzNum']-1] = '✓'
            local_storage.setItem(itemKey='history',itemValue=hist1,key='win')
        else:
            st.title("Nope.", text_alignment='center')
            hist1 = local_storage.getItem('history')
            hist1[st.session_state['puzzNum']-1] = '✘'
            local_storage.setItem(itemKey='history',itemValue=hist1,key='lose')
        st.subheader(st.session_state.key.at['Name','Value'],text_alignment = 'center')
        
        if expH:
            st.image(st.session_state.imageKeyH)
        else:
            st.image(st.session_state.imageKey)
        st.text(st.session_state.key.at['Info','Value'], text_alignment = 'justify')

st.caption('© 2026 Kelsey Sakimoto', text_alignment='right')
