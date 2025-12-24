"""
Database utility functions for RXNRECer
Supports MySQL, PostgreSQL, and SQLite databases
"""
import os
import logging
from typing import Optional, Dict, Any, List
from contextlib import contextmanager
import sqlite3

logger = logging.getLogger(__name__)

# Try to import database drivers
try:
    import pymysql
    MYSQL_AVAILABLE = True
except ImportError:
    MYSQL_AVAILABLE = False
    logger.warning("PyMySQL not available. MySQL support disabled.")

try:
    import psycopg2
    from psycopg2 import pool
    POSTGRESQL_AVAILABLE = True
except ImportError:
    POSTGRESQL_AVAILABLE = False
    logger.warning("psycopg2 not available. PostgreSQL support disabled.")


class DatabaseConnection:
    """Database connection manager supporting multiple database types"""
    
    def __init__(self, db_type: str = 'sqlite', **kwargs):
        """
        Initialize database connection
        
        Args:
            db_type: Database type ('mysql', 'postgresql', 'sqlite')
            **kwargs: Database connection parameters
                For MySQL: host, port, user, password, database
                For PostgreSQL: host, port, user, password, database
                For SQLite: database (path to .db file)
        """
        self.db_type = db_type.lower()
        self.connection_params = kwargs
        self._connection = None
        self._connection_pool = None
        
        if self.db_type == 'mysql':
            if not MYSQL_AVAILABLE:
                raise ImportError("PyMySQL is required for MySQL support. Install it with: pip install pymysql")
            self._init_mysql()
        elif self.db_type == 'postgresql' or self.db_type == 'postgres':
            if not POSTGRESQL_AVAILABLE:
                raise ImportError("psycopg2 is required for PostgreSQL support. Install it with: pip install psycopg2-binary")
            self._init_postgresql()
        elif self.db_type == 'sqlite':
            self._init_sqlite()
        else:
            raise ValueError(f"Unsupported database type: {db_type}. Supported types: mysql, postgresql, sqlite")
    
    def _init_mysql(self):
        """Initialize MySQL connection"""
        required_params = ['host', 'user', 'password', 'database']
        for param in required_params:
            if param not in self.connection_params:
                raise ValueError(f"Missing required parameter for MySQL: {param}")
        
        self.connection_params.setdefault('port', 3306)
        self.connection_params.setdefault('charset', 'utf8mb4')
        self.connection_params.setdefault('autocommit', False)
    
    def _init_postgresql(self):
        """Initialize PostgreSQL connection"""
        required_params = ['host', 'user', 'password', 'database']
        for param in required_params:
            if param not in self.connection_params:
                raise ValueError(f"Missing required parameter for PostgreSQL: {param}")
        
        self.connection_params.setdefault('port', 5432)
    
    def _init_sqlite(self):
        """Initialize SQLite connection"""
        if 'database' not in self.connection_params:
            raise ValueError("Missing required parameter for SQLite: database (path to .db file)")
    
    def connect(self):
        """Establish database connection"""
        if self._connection:
            return self._connection
        
        try:
            if self.db_type == 'mysql':
                self._connection = pymysql.connect(**self.connection_params)
            elif self.db_type == 'postgresql' or self.db_type == 'postgres':
                self._connection = psycopg2.connect(**self.connection_params)
            elif self.db_type == 'sqlite':
                self._connection = sqlite3.connect(
                    self.connection_params['database'],
                    check_same_thread=False
                )
                self._connection.row_factory = sqlite3.Row
            
            logger.info(f"Successfully connected to {self.db_type} database")
            return self._connection
        except Exception as e:
            logger.error(f"Failed to connect to database: {e}")
            raise
    
    def disconnect(self):
        """Close database connection"""
        if self._connection:
            try:
                if self.db_type == 'sqlite':
                    self._connection.close()
                else:
                    self._connection.close()
                self._connection = None
                logger.info("Database connection closed")
            except Exception as e:
                logger.error(f"Error closing database connection: {e}")
    
    @contextmanager
    def get_connection(self):
        """Context manager for database connections"""
        conn = self.connect()
        try:
            yield conn
        finally:
            if self.db_type != 'sqlite':  # SQLite connections are kept open
                pass  # Don't close for connection pooling
    
    def execute_query(self, query: str, params: Optional[tuple] = None) -> List[Dict[str, Any]]:
        """
        Execute a SELECT query and return results
        
        Args:
            query: SQL query string
            params: Query parameters (for parameterized queries)
            
        Returns:
            List of dictionaries representing rows
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            try:
                if params:
                    cursor.execute(query, params)
                else:
                    cursor.execute(query)
                
                # Fetch results
                if self.db_type == 'sqlite':
                    rows = cursor.fetchall()
                    # Convert sqlite3.Row to dict
                    result = [dict(row) for row in rows]
                else:
                    # For MySQL and PostgreSQL
                    columns = [desc[0] for desc in cursor.description] if cursor.description else []
                    rows = cursor.fetchall()
                    result = [dict(zip(columns, row)) for row in rows]
                
                return result
            except Exception as e:
                logger.error(f"Error executing query: {e}")
                raise
            finally:
                cursor.close()
    
    def execute_update(self, query: str, params: Optional[tuple] = None) -> int:
        """
        Execute an INSERT, UPDATE, or DELETE query
        
        Args:
            query: SQL query string
            params: Query parameters (for parameterized queries)
            
        Returns:
            Number of affected rows
        """
        with self.get_connection() as conn:
            cursor = conn.cursor()
            try:
                if params:
                    cursor.execute(query, params)
                else:
                    cursor.execute(query)
                
                if self.db_type != 'sqlite':
                    conn.commit()
                else:
                    conn.commit()
                
                return cursor.rowcount
            except Exception as e:
                if self.db_type != 'sqlite':
                    conn.rollback()
                else:
                    conn.rollback()
                logger.error(f"Error executing update: {e}")
                raise
            finally:
                cursor.close()
    
    def __enter__(self):
        """Context manager entry"""
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.disconnect()


def get_db_connection_from_env() -> DatabaseConnection:
    """
    Create database connection from environment variables
    
    Environment variables:
        DB_TYPE: Database type (mysql, postgresql, sqlite)
        DB_HOST: Database host (for MySQL/PostgreSQL)
        DB_PORT: Database port (for MySQL/PostgreSQL)
        DB_USER: Database user (for MySQL/PostgreSQL)
        DB_PASSWORD: Database password (for MySQL/PostgreSQL)
        DB_NAME: Database name (for MySQL/PostgreSQL/SQLite)
        DB_PATH: Database file path (for SQLite, alternative to DB_NAME)
    """
    db_type = os.environ.get('DB_TYPE', 'sqlite').lower()
    
    if db_type == 'mysql':
        return DatabaseConnection(
            db_type='mysql',
            host=os.environ.get('DB_HOST', 'localhost'),
            port=int(os.environ.get('DB_PORT', 3306)),
            user=os.environ.get('DB_USER'),
            password=os.environ.get('DB_PASSWORD'),
            database=os.environ.get('DB_NAME')
        )
    elif db_type in ('postgresql', 'postgres'):
        return DatabaseConnection(
            db_type='postgresql',
            host=os.environ.get('DB_HOST', 'localhost'),
            port=int(os.environ.get('DB_PORT', 5432)),
            user=os.environ.get('DB_USER'),
            password=os.environ.get('DB_PASSWORD'),
            database=os.environ.get('DB_NAME')
        )
    else:  # sqlite
        db_path = os.environ.get('DB_PATH') or os.environ.get('DB_NAME', 'rxnrecer_config.db')
        return DatabaseConnection(
            db_type='sqlite',
            database=db_path
        )

