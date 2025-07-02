document.getElementById('calculate-btn').addEventListener('click', function() {
    const segmentUrl = document.getElementById('segment-url').value;
    const profileUrl = document.getElementById('profile-url').value;
    
    if (!segmentUrl || !profileUrl) {
        alert('Please enter both URLs');
        return;
    }
    
    // Here you would add your logic to:
    // 1. Extract segment ID and athlete ID from URLs
    // 2. Make API calls to Strava (you'll need authentication)
    // 3. Process the data
    // 4. Display results
    
    // For now, just a placeholder
    const resultsDiv = document.getElementById('results');
    resultsDiv.innerHTML = `
        <h3>Results</h3>
        <p>Segment: ${segmentUrl}</p>
        <p>Athlete: ${profileUrl}</p>
        <p>This is where your performance data would appear.</p>
    `;
});
