let currentOutputDir = null; 
var jmolApplet0;

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
    $("#viewer-container").html(Jmol.getAppletHtml("jmolApplet0", Info));
    setStatus("System initialized. Awaiting directories...");
    try {
        let res = await fetch('/api/defaults');
        let data = await res.json();
        if (data.nodes) document.getElementById('nodes-path').value = data.nodes;
        if (data.linkers) document.getElementById('linkers-path').value = data.linkers;
        if (data.out) document.getElementById('out-path').value = data.out;
        if (data.nodes && data.linkers && data.out) loadDirectories();
    } catch (e) {
        setStatus("Could not fetch default paths.", true);
    }
};

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
            loadMolecule(fileUrl);
        };
        historyList.prepend(item); 
    }
}

async function loadDirectories() {
    let nodes = document.getElementById('nodes-path').value.trim();
    let linkers = document.getElementById('linkers-path').value.trim();
    let out = document.getElementById('out-path').value.trim();

    if (!nodes || !linkers || !out) {
        setStatus("Please fill in all three directory paths.", true);
        return;
    }

    setStatus("Loading directories...");
    document.getElementById('btn-load-dirs').disabled = true;

    try {
        let response = await fetch('/api/config', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ nodes: nodes, linkers: linkers, out: out })
        });

        let data = await response.json();
        if (!response.ok) throw new Error(data.error || "Failed to load directories");

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
            setStatus(`Loaded ${data.history.length} previous assemblies.`);
        }
        document.getElementById('btn-build').disabled = false;
        setStatus(data.message);
    } catch (error) {
        setStatus("Error: " + error.message, true);
    } finally {
        document.getElementById('btn-load-dirs').disabled = false;
    }
}

function populateSelect(elementId, items) {
    let select = document.getElementById(elementId);
    select.innerHTML = '';
    if (!items || items.length === 0) {
        select.innerHTML = '<option value="">None found</option>';
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
            body: JSON.stringify({ node: node, linker: linker, shape: shape })
        });

        let data = await response.json();
        if (!response.ok) throw new Error(data.error || "Unknown build error");

        setStatus("Build successful! Rendering...");
        currentOutputDir = data.output_dir;
        document.getElementById('btn-optimize').disabled = false; 
        
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

    // Disable buttons
    btnBuild.disabled = true;
    btnOpt.disabled = true;

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
                    body: JSON.stringify({ target_dir: currentOutputDir, method: method })
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
                btnOpt.innerHTML = originalText;
            }
        });
    });
}

function loadMolecule(fileUrl) {
    try {
        let cacheBusterUrl = fileUrl + "?t=" + new Date().getTime();
        let script = `
            load "${cacheBusterUrl}";
            wireframe 0.15;
            spacefill 20%;
            zoom center;
        `;
        Jmol.script(jmolApplet0, script);
        setStatus("Rendering complete.");
    } catch (error) {
        setStatus("Failed to load molecule viewer: " + error.message, true);
    }
}
