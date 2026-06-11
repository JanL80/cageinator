let currentOutputDir = null; 
var jmolApplet0;
let advancedOpen = false;

function getSessionId() {
    let sid = localStorage.getItem('cageinator_session');
    if (!sid) {
        sid = Math.random().toString(36).substr(2, 16);
        localStorage.setItem('cageinator_session', sid);
    }
    return sid;
}
const SESSION_ID = getSessionId();

var Info = {
    width: "100%",
    height: "100%",
    debug: false,
    color: "#1e1e1e",           
    addSelectionOptions: false,
    use: "HTML5",
    j2sPath: "/static/jsmol/j2s", 
    disableJ2SLoadMonitor: true,
    disableInitialConsole: true
};

window.onload = async function() {
    Jmol.setDocument(0);     
    jmolApplet0 = Jmol.getApplet("jmolApplet0", Info); 
    $("#viewer-container").html(jmolApplet0._code); 
    setStatus("System initialized. Restoring session...");
    document.getElementById('upload-node-file').addEventListener('change', (e) => uploadFile(e, 'nodes'));
    document.getElementById('upload-linker-file').addEventListener('change', (e) => uploadFile(e, 'linkers'));
    restoreSession();
};

async function restoreSession() {
    try {
        let response = await fetch('/api/sync_session', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: SESSION_ID })
        });

        let data = await response.json();
        if (!response.ok) throw new Error(data.error || "Failed to sync session");

        populateSelect('linker-select', data.linkers);
        populateSelect('node-select', data.nodes);
        populateSelect('shape-select', data.shapes);
        
        document.getElementById('history-list').innerHTML = ''; 
        if (data.history && data.history.length > 0) {
            for (let i = data.history.length - 1; i >= 0; i--) {
                let item = data.history[i];
                let tag = item.tag || "Restored Assembly";
                addHistoryItem(item.name, item.url, tag);
            }
        }
        
        if(data.linkers.length > 0 && data.nodes.length > 0) {
            document.getElementById('btn-build').disabled = false;
        }
        
        setStatus(data.message);
    } catch (error) {
        setStatus("Session Warning: " + error.message, true);
    }
}

const svgMenu = `<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><line x1="3" y1="12" x2="21" y2="12"></line><line x1="3" y1="6" x2="21" y2="6"></line><line x1="3" y1="18" x2="21" y2="18"></line></svg>`;
const svgClose = `<svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="#c00505" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><line x1="18" y1="6" x2="6" y2="18"></line><line x1="6" y1="6" x2="18" y2="18"></line></svg>`;

function toggleAdvancedMenu() {
    advancedOpen = !advancedOpen;
    let overlay = document.getElementById('advanced-overlay');
    let btn = document.getElementById('advanced-menu');
    
    if (advancedOpen) {
        overlay.style.display = 'flex';
        btn.innerHTML = svgClose;
    } else {
        overlay.style.display = 'none';
        btn.innerHTML = svgMenu;
    }
}

async function uploadFile(event, type) {
    let file = event.target.files[0];
    if (!file) return;

    let formData = new FormData();
    formData.append('file', file);
    formData.append('type', type);
    formData.append('session_id', SESSION_ID);

    setStatus(`Uploading ${file.name}...`);
    try {
        let response = await fetch('/api/upload', {
            method: 'POST',
            body: formData
        });
        let data = await response.json();
        if (!response.ok) throw new Error(data.error || "Upload failed");

        populateSelect('linker-select', data.linkers);
        populateSelect('node-select', data.nodes);
        populateSelect('shape-select', data.shapes);
        
        document.getElementById('btn-build').disabled = false;
        setStatus(data.message);
    } catch (err) {
        setStatus("Error: " + err.message, true);
    }
    event.target.value = ''; 
}

function setStatus(msg, isError = false) {
    let logBox = document.getElementById('status-box');
    let entry = document.createElement('div');
    entry.className = 'log-entry' + (isError ? ' error' : '');
    
    let now = new Date();
    let timeString = now.toLocaleTimeString([], {hour: '2-digit', minute:'2-digit', second:'2-digit'});
    
    entry.innerHTML = `<span class="log-time">[${timeString}]</span> ${msg}`;
    logBox.appendChild(entry);
    logBox.scrollTop = logBox.scrollHeight;
}

function addHistoryItem(dirName, fileUrl, type) {
    let historyList = document.getElementById('history-list');
    let safeId = 'hist-' + dirName.replace(/[^a-zA-Z0-9_-]/g, '-');
    let existingItem = document.getElementById(safeId);

    if (existingItem) {
        existingItem.querySelector('.history-item-subtitle').innerText = type;
        existingItem.onclick = () => {
            setStatus(`Reloading past assembly: ${dirName}`);
            currentOutputDir = dirName; 
            document.getElementById('btn-optimize').disabled = false;
            document.getElementById('btn-download').disabled = false;
            loadMolecule(fileUrl);
        };
        historyList.prepend(existingItem);
    } else {
        let item = document.createElement('div');
        item.className = 'history-item';
        item.id = safeId;
        item.innerHTML = `
            <div class="history-item-title" title="${dirName}">${dirName}</div>
            <div class="history-item-subtitle">${type}</div>
        `;
        item.onclick = () => {
            setStatus(`Reloading past assembly: ${dirName}`);
            currentOutputDir = dirName; 
            document.getElementById('btn-optimize').disabled = false;
            document.getElementById('btn-download').disabled = false;
            loadMolecule(fileUrl);
        };
        historyList.prepend(item); 
    }
}

function populateSelect(elementId, items) {
    let select = document.getElementById(elementId);
    select.innerHTML = '';
    if (!items || items.length === 0) {
        select.innerHTML = '<option value="">None available</option>';
        return;
    }
    items.forEach(item => {
        let option = document.createElement('option');
        option.value = item;
        option.innerText = item;
        select.appendChild(option);
    });
}

async function buildCage() {
    let node = document.getElementById('node-select').value;
    let linker = document.getElementById('linker-select').value;
    let shape = document.getElementById('shape-select').value;

    if (!node || !linker || !shape) {
        setStatus("Please select a Node, Linker, and Shape.", true);
        return;
    }

    setStatus(`Assembling ${node} + ${linker} into ${shape}...`);
    document.getElementById('btn-build').disabled = true;

    try {
        let response = await fetch('/api/build', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: SESSION_ID, node: node, linker: linker, shape: shape })
        });

        let data = await response.json();
        if (!response.ok) throw new Error(data.error || "Unknown build error");

        setStatus("Build successful! Rendering...");
        currentOutputDir = data.output_dir;
        
        document.getElementById('btn-optimize').disabled = false; 
        document.getElementById('btn-download').disabled = false;
        
        addHistoryItem(data.output_dir, data.file_url, "Assembled");
        loadMolecule(data.file_url);
    } catch (error) {
        setStatus("Build Failed: " + error.message, true);
    } finally {
        document.getElementById('btn-build').disabled = false;
    }
}

function optimizeCage() {
    if (!currentOutputDir) return;

    let method = document.getElementById('opt-method').value;
    setStatus(`Optimizing with ${method.toUpperCase()}... This may take a while.`);
    
    let btnBuild = document.getElementById('btn-build');
    let btnOpt = document.getElementById('btn-optimize');
    let btnDl = document.getElementById('btn-download');

    // Disable buttons
    btnBuild.disabled = true;
    btnOpt.disabled = true;
    btnDl.disabled = true;

    // Inject spinner directly into the button text
    let originalText = btnOpt.innerHTML;
    btnOpt.innerHTML = `<span style="display: flex; align-items: center; justify-content: center; gap: 8px;">Optimizing... <span class="spinner"></span></span>`;

    // Force browser to paint the button before locking the thread
    requestAnimationFrame(() => {
        requestAnimationFrame(async () => {
            try {
                let response = await fetch('/api/optimize', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ session_id: SESSION_ID, target_dir: currentOutputDir, method: method })
                });

                let data = await response.json();
                if (!response.ok) throw new Error(data.error || "Unknown optimization error");

                setStatus(`Optimization successful! Reloading structure...`);
                addHistoryItem(currentOutputDir, data.file_url, `Optimized (${method.toUpperCase()})`);
                loadMolecule(data.file_url);
                
            } catch (error) {
                setStatus("Optimization Failed: " + error.message, true);
            } finally {
                // Restore buttons
                btnBuild.disabled = false;
                btnOpt.disabled = false;
                btnDl.disabled = false;
                btnOpt.innerHTML = originalText;
            }
        });
    });
}

function triggerDownload(url) {
    let iframe = document.createElement('iframe');
    iframe.style.display = 'none';
    iframe.src = url;
    document.body.appendChild(iframe);
    
    setTimeout(() => {
        if (document.body.contains(iframe)) {
            document.body.removeChild(iframe);
        }
    }, 10000);
}

function downloadCage() {
    if (!currentOutputDir) return;
    setStatus(`Preparing ZIP file for ${currentOutputDir}...`);
    
    let url = `/api/download?session_id=${SESSION_ID}&target_dir=${encodeURIComponent(currentOutputDir)}`;
    triggerDownload(url);
}

function downloadWorkspace() {
    setStatus(`Preparing full workspace ZIP...`);
    
    let url = `/api/download_workspace?session_id=${SESSION_ID}`;
    triggerDownload(url);
}

async function clearSession() {
    if(!confirm("Are you sure you want to clear your session? This will wipe all uploaded files and built cages.")) return;
    
    try {
        let res = await fetch('/api/clear_session', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ session_id: SESSION_ID })
        });
        
        if(!res.ok) throw new Error("Failed to clear backend session.");
        
        document.getElementById('history-list').innerHTML = '';
        document.getElementById('status-box').innerHTML = '';
        populateSelect('linker-select', []);
        populateSelect('node-select', []);
        
        if (jmolApplet0) Jmol.script(jmolApplet0, "zap;");
        
        currentOutputDir = null;
        document.getElementById('btn-optimize').disabled = true;
        document.getElementById('btn-download').disabled = true;
        document.getElementById('btn-build').disabled = true;
        
        setStatus("Session wiped clean. Ready for new uploads.");
    } catch (e) {
        setStatus("Error clearing session: " + e.message, true);
    }
}

async function loadMolecule(fileUrl) {
    try {
        let cacheBusterUrl = fileUrl + "?t=" + new Date().getTime();
        let response = await fetch(cacheBusterUrl);
        
        if (!response.ok) {
            throw new Error("Could not retrieve structure data from the server.");
        }
        
        let molData = await response.text();
        
        let script = `
            zap;
            load DATA "inline_model"
${molData}
            END "inline_model";
            wireframe 0.15;
            spacefill 20%;
            zoom center;
            refresh;
        `;
        
        Jmol.script(jmolApplet0, script);
        setStatus("Rendering complete.");
    } catch (error) {
        setStatus("Failed to load molecule viewer: " + error.message, true);
    }
}
