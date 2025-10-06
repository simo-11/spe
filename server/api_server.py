# -*- coding: utf-8 -*-
"""
Created on Mon Oct  6 20:22:26 2025

@author: simon and copilot
"""
from fastapi import FastAPI
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from analyze import analyze

app = FastAPI()

class AnalysisRequest(BaseModel):
    param1: float
    param2: float
    param3: float
    param4: float
    param5: float

@app.post("/analyze")
async def analyze_endpoint(request: AnalysisRequest):
    params = request.dict()
    try:
        return await analyze(params)
    except ValueError as e:
        return JSONResponse(status_code=422, content={"error": str(e)})
    
if __name__ == "__main__":
    import uvicorn
    uvicorn.run("api_server:app", host="127.0.0.1", port=8000, reload=False)
